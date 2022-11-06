library(parallel)
library(tidyverse)
library(simulateGP)
library(pROC)
library(pbapply)
library(lavaan)

prop_overlap = function(b_disc, b_rep, se_disc, se_rep, alpha)
{
  p_sign <- pnorm(-abs(b_disc) / se_disc) * pnorm(-abs(b_disc) / se_rep) + ((1 - pnorm(-abs(b_disc) / se_disc)) * (1 - pnorm(-abs(b_disc) / se_rep)))
  p_sig <- pnorm(-abs(b_disc) / se_rep + qnorm(alpha / 2)) + (1 - pnorm(-abs(b_disc) / se_rep - qnorm(alpha / 2)))
  p_rep <- pnorm(abs(b_rep)/se_rep, lower.tail=FALSE)
  res <- dplyr::tibble(
    nsnp=length(b_disc),
    metric=c("Sign", "Sign", "P-value", "P-value"),
    datum=c("Expected", "Observed", "Expected", "Observed"),
    value=c(sum(p_sign, na.rm=TRUE), sum(sign(b_disc) == sign(b_rep)), sum(p_sig, na.rm=TRUE), sum(p_rep < alpha, na.rm=TRUE))
  )
  return(list(res=res, variants=dplyr::tibble(sig=p_sig, sign=p_sign, fdisc=b_disc^2/se_disc^2, frep=b_rep^2/se_rep^2)))
}


organise_ldobj <- function(dir1, dir2, region)
{
  ld1 <- readRDS(file.path(dir1, paste0("ldobj_", region, ".rds")))
  ld2 <- readRDS(file.path(dir2, paste0("ldobj_", region, ".rds")))
  x1 <- unique(ld1$map$snp)
  x1 <- x1[x1 %in% ld1$map$snp]
  x1 <- x1[x1 %in% ld2$map$snp]
  i1 <- which(ld1$map$snp %in% x1)
  ld1$map <- ld1$map[i1,]
  ld1$ld <- ld1$ld[i1, i1]
  i2 <- which(ld2$map$snp %in% x1)
  ld2$map <- ld2$map[i2,]
  ld2$ld <- ld2$ld[i2, i2]
  stopifnot(all(ld1$map$snp == ld2$map$snp))
  return(list(ld1=ld1, ld2=ld2))
}


subset_ldobj <- function(ldobj, snps)
{
  i <- which(ldobj$map$snp %in% snps)
  ldobj$map <- ldobj$map[i,]
  ldobj$ld <- ldobj$ld[i,i]
  return(ldobj)
}


# 1. Generate a causal variant for pop1
# 2. keep the region around the causal variant
# 3. generate LD
# 4. 

simss <- function(ld1, ld2, nid1, nid2, type, hsq1, window, bxy1, bxy2)
{
  x1p <- ld1$map %>%
    generate_gwas_params(h2=hsq1, S=0, Pi=1/nrow(.))
  pos1 <- x1p$pos[x1p$beta != 0] - window/2
  pos2 <- pos1 + window
  
  x1p <- subset(x1p, pos >= pos1 & pos <= pos2)  
  ld1 <- subset_ldobj(ld1, x1p$snp)
  ld2 <- subset_ldobj(ld2, x1p$snp)
  
  x1 <- x1p %>%
    generate_gwas_ss(nid=nid1, ldobj=ld1)

  y1p <- x1p
  y1p$beta <- y1p$beta * bxy1
  y1 <- y1p %>%
    generate_gwas_ss(nid=nid1, ldobj=ld1)
  
  x2p <- ld2$map %>% generate_gwas_params(h2=hsq1, S=0, Pi=1/nrow(.))
  if(type=="shared")
  {
    x2p$beta <- x1p$beta
  } else if(type == "1") {
    x2p$beta <- 0
  }

  x2 <- x2p %>%
    generate_gwas_ss(nid=nid2, ldobj=ld2)
  
  y2p <- x2p
  y2p$beta <- y2p$beta * bxy2
  y2 <- y2p %>%
    generate_gwas_ss(nid=nid2, ldobj=ld2)

  return(list(x1=x1, x2=x2, y1=y1, y2=y2))
}


plot_simss <- function(s)
{
  raw <- which.max(s[[1]]$fval)
  maxz <- which.max(s[[1]]$fval + s[[2]]$fval)
  i <- tibble(pos=c(s[[1]]$pos[raw], s[[1]]$pos[maxz]), method=c("raw", "maxz"))
  rbind(
    tibble(snp=s[[1]]$snp, pos=s[[1]]$pos, pval=-log10(s[[1]]$pval), causal=s[[1]]$beta != 0, what="pop1"),
    tibble(snp=s[[1]]$snp, pos=s[[1]]$pos, pval=-log10(s[[2]]$pval), causal=s[[2]]$beta != 0, what="pop2")
  ) %>%
    ggplot(., aes(x=pos, y=pval)) +
    geom_point(aes(colour=causal, size=causal, shape=what)) +
    facet_grid(what ~ .) +
    geom_hline(yintercept=0) +
    geom_vline(data=i, aes(xintercept=pos, linetype=method))
}

runsem <- function(model, data, modname)
{
  mod <- lavaan::sem(model, data=data)
  mod <- lavaan::summary(mod, fit.measures=TRUE)
  o <- tibble::tibble(
                      Methods=modname,
                      pop=1:2,
                      nsnp=nrow(data),
                      bivhat=mod$PE$est[1:2],
                      se=mod$PE$se[1:2],
                      pval=mod$PE$pvalue[1:2],
                      aic=mod$FIT['aic']
                    ) %>%  dplyr::mutate(pop = as.character(pop))
  return(o)
}

perform_basic_sem = function(harmonised_dat)
{
  d <- harmonised_dat %>%
   dplyr::mutate(r1 = y1/x1) %>%
   dplyr::mutate(r2 = y2/x2) %>%
   dplyr::mutate(w1 = sqrt(x1^2 / yse1^2)) %>%
   dplyr::mutate(w2 = sqrt(x2^2 / yse2^2)) %>%
   dplyr::mutate(o1 = r1 * w1) %>%
   dplyr::mutate(o2 = r2 * w2)

  out <- list()
  out$ivw1 <- TwoSampleMR::mr_ivw(d$x1, d$y1, d$xse1, d$yse1) %>%
  {tibble::tibble(Methods="IVW", pop="1", nsnp=nrow(d), bivhat=.$b, se=.$se, pval=.$pval)}
  out$ivw2 <- TwoSampleMR::mr_ivw(d$x2, d$y2, d$xse2, d$yse2) %>%
  {tibble::tibble(Methods="IVW", pop="2", nsnp=nrow(d), bivhat=.$b, se=.$se, pval=.$pval)}
  out$rm1 <- summary(lm(o1 ~ -1 + w1, data=d)) %>%
  {tibble::tibble(Methods="RadialIVW", pop="1", nsnp=nrow(d), bivhat=.$coef[1,1], se=.$coef[1,2], pval=.$coef[1,4])}
  out$rm2 <- summary(lm(o2 ~ -1 + w2, data=d)) %>%
  {tibble::tibble(Methods="RadialIVW", pop="2", nsnp=nrow(d), bivhat=.$coef[1,1], se=.$coef[1,2], pval=.$coef[1,4])}

  out$semA <- runsem('
  y1 ~ biv*x1
  y2 ~ biv*x2
  ', d, "UnweightedSEMa")[1, ] %>%
  dplyr::mutate(pop=replace(pop, pop==1, "1=2"))

  out$semB <- runsem('
  y1 ~ biv_1*x1
  y2 ~ biv_2*x2
  ', d, "UnweightedSEMb")

  out$modC <- runsem('
  o1 ~ biv*w1
  o2 ~ biv*w2
  ', d, "RadialSEMa")[1, ] %>%
  dplyr::mutate(pop=replace(pop, pop==1, "1=2"))

  out$modD <- runsem('
  o1 ~ biv_1*w1
  o2 ~ biv_2*w2
  ', d, "RadialSEMb")

  return(dplyr::bind_rows(out))
}


sim <- function(lddir1, lddir2, region, nid1, nid2, nsim=100, pshared, pdistinct, p1, hsq1=0.4, window=250000, winnerscurse=TRUE, removecv=TRUE, bxy1, bxy2, sim=1, simid=1, mc.cores=1)
{
  type <- sample(c("shared", "distinct", "1"), nsim, replace=T, prob=c(pshared, pdistinct, p1))
  ld <- organise_ldobj(lddir1, lddir2, region)
  hsq <- 1:nsim / sum(1:nsim) * hsq1
  l <- pblapply(1:nsim, function(i)
  {
    x <- simss(ld[[1]], ld[[2]], nid1, nid2, type[i], hsq[i], window, bxy1, bxy2)
    x1 <- x$x1
    if(!any(x1$pval < 5e-8))
    {
      return(NULL)
    }
    x2 <- x$x2
    if(removecv)
    {
      cv <- which(x1$beta != 0)
      x1 <- x1[-cv,]
      x2 <- x2[-cv,]
    }
    # Get the top hit for 1 and check for overlap in 2
    x1raw <- subset(x1, fval==max(fval, na.rm=T))[1,]
    x2raw <- subset(x2, snp == x1raw$snp)
    if(!winnerscurse)
    {
      x1raw$bhat <- rnorm(nrow(x1raw), x1raw$beta, x1raw$se)
    }
    oraw <- prop_overlap(x1raw$bhat, x2raw$bhat, x1raw$se, x2raw$se, 0.05/nsim)

    fvalsum <- x1$fval + x2$fval
    maxz <- which.max(fvalsum)[1]
    x1maxz <- x1[maxz, ]
    x2maxz <- x2[maxz, ]
    if(!winnerscurse)
    {
      x1maxz$bhat <- rnorm(nrow(x1maxz), x1maxz$beta, x1maxz$se)
    }
    omaxz <- prop_overlap(x1maxz$bhat, x2maxz$bhat, x2maxz$se, x2maxz$se, 0.05/nsim)

    y1raw <- subset(x$y1, snp == x1raw$snp)
    y2raw <- subset(x$y2, snp == x1raw$snp)
    y1maxz <- subset(x$y1, snp == x1maxz$snp)
    y2maxz <- subset(x$y2, snp == x1maxz$snp)

    datraw <- tibble(SNP=i, x1=x1raw$bhat, y1=y1raw$bhat, xse1=x1raw$se, yse1=y1raw$se, x2=x2raw$bhat, y2=y2raw$bhat, xse2=x2raw$se, yse2=y2raw$se)
    datmaxz <- tibble(SNP=i, x1=x1maxz$bhat, y1=y1maxz$bhat, xse1=x1maxz$se, yse1=y1maxz$se, x2=x2maxz$bhat, y2=y2maxz$bhat, xse2=x2maxz$se, yse2=y2maxz$se)

    return(list(raw=oraw, maxz=omaxz, th1=x1raw$snp, th2=x1maxz$snp, datraw=datraw, datmaxz=datmaxz, x2pval=x2raw$pval))
  }, cl=mc.cores)
  
  sig <- which(sapply(l, function(x) !is.null(x)))
  l[sapply(l, is.null)] <- NULL
  prop_same <- sapply(l, function(x) x$th1 == x$th2) %>% sum %>% {./length(l)}
  o1 <- list(res=l[[1]]$raw$res)
  o1$variants <- lapply(l, function(x)
  {
    x$raw$variants
  }) %>% bind_rows()
  o1$res$value <- sapply(l, function(x) x$raw$res$value) %>% rowSums()
  
  o1$variants$snp <- sig
  o1$variants$truth <- type[sig] %in% c("shared", "distinct")
  o1$variants$obs_sig <- sapply(l, function(x) x$raw$res$value[4])
  o1$variants$obs_sign <- sapply(l, function(x) x$raw$res$value[2])
  o1$variants$correct_sig <- o1$variants$truth == o1$variants$obs_sig
  o1$variants$pred <- o1$variants$sig
  o1$variants$pred[!o1$variants$obs_sig] <- 1 - o1$variants$pred[!o1$variants$obs_sig]


  o2 <- list(res=l[[1]]$maxz$res)
  o2$variants <- lapply(l, function(x)
  {
    x$maxz$variants
  }) %>% bind_rows()
  o2$res$value <- sapply(l, function(x) x$maxz$res$value) %>% rowSums()
  
  o2$variants$snp <- sig
  o2$variants$truth <- type[sig] %in% c("shared", "distinct")
  o2$variants$obs_sig <- sapply(l, function(x) x$maxz$res$value[4])
  o2$variants$obs_sign <- sapply(l, function(x) x$maxz$res$value[2])
  o2$variants$correct_sig <- o2$variants$truth == o2$variants$obs_sig
  o2$variants$pred <- o2$variants$sig
  o2$variants$pred[!o2$variants$obs_sig] <- 1 - o2$variants$pred[!o2$variants$obs_sig]
  
  # truth = whether a variant actually has an effect in both pops
  # pred = prob(sig) (for overlaps) or 1-prob(sig) (for non-overlaps)
  # auc = extent to which pred can predict truth
  o1$auc <- tryCatch(as.numeric(suppressMessages(auc(o1$variants$truth, o1$variants$pred))), error=function(e) return(NA))
  o2$auc <- tryCatch(as.numeric(suppressMessages(auc(o2$variants$truth, o2$variants$pred))), error=function(e) return(NA))
  
  perf1 <- o1$variants %>% group_by(truth) %>% summarise(pow=sum(frep > qf(0.05/nsim, 1, nid2, low=F))/n(), meanf=mean(frep))
  perf2 <- o2$variants %>% group_by(truth) %>% summarise(pow=sum(frep > qf(0.05/nsim, 1, nid2, low=F))/n(), meanf=mean(frep))
  
  obs_rep1 = subset(o1$res, metric == "P-value" & datum == "Observed")$value[1]
  obs_rep2 = subset(o2$res, metric == "P-value" & datum == "Observed")$value[1]
  exp_rep1 = subset(o1$res, metric == "P-value" & datum == "Expected")$value[1]
  exp_rep2 = subset(o2$res, metric == "P-value" & datum == "Expected")$value[1]
  binom1 = binom.test(x=obs_rep1, n=length(sig), p=exp_rep1/length(sig))
  binom2 = binom.test(x=obs_rep2, n=length(sig), p=exp_rep2/length(sig))
  
  # mr
  datraw <- lapply(l, function(x) x$datraw) %>% bind_rows()
  datmaxz <- lapply(l, function(x) x$datmaxz) %>% bind_rows()
  mrres <- bind_rows(
    perform_basic_sem(datraw) %>% mutate(instruments="raw"),
    perform_basic_sem(datmaxz) %>% mutate(instruments="maxz")
  ) %>% mutate(simid=simid)

  res <- tibble(
    lddir1=lddir1,
    lddir2=lddir2,
    region=region,
    nid1=nid1,
    nid2=nid2,
    pshared=pshared,
    pdistinct=pdistinct,
    p1=p1,
    hsq1=hsq1,
    window=window,
    winnerscurse=winnerscurse,
    nsig=length(sig),
    nsig2=sum(sapply(l, function(x) x$x2pval) < (0.05/nsig)),
    method=c("raw", "maxz"),
    auc=c(o1$auc, o2$auc),
    sim=sim,
    simid=simid,
    correct_sig=
      c(sum(o1$variants$correct_sig)/length(o1$variants$correct_sig),
      sum(o2$variants$correct_sig)/length(o2$variants$correct_sig)),
    obs_rep=c(obs_rep1,obs_rep2),
    exp_rep=c(exp_rep1, exp_rep2),
    binom_p = c(binom1$p.value, binom2$p.value),
    fdr=c(subset(perf1, !truth)$pow,subset(perf2, !truth)$pow),
    pow=c(subset(perf1, truth)$pow,subset(perf2, truth)$pow),
    prop_same=prop_same
  )
  
  return(list(instruments=res, mr=mrres))
}
