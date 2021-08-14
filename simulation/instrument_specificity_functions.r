library(parallel)
library(tidyverse)
library(simulateGP)
library(pROC)
library(pbapply)

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


simss <- function(ld1, ld2, nid1, nid2, type, hsq1, window)
{
  x1 <- ld1$map %>%
    generate_gwas_params(h2=hsq1, S=0, Pi=1/nrow(.))
  pos1 <- x1$pos[x1$beta != 0] - window/2
  pos2 <- pos1 + window
  
  x1 <- subset(x1, pos >= pos1 & pos <= pos2)  
  ld1 <- subset_ldobj(ld1, x1$snp)
  ld2 <- subset_ldobj(ld2, x1$snp)
  
  x1 <- x1 %>%
    generate_gwas_ss(nid=nid1, ldobj=ld1)
  
  x2 <- ld2$map %>% generate_gwas_params(h2=hsq1, S=0, Pi=1/nrow(.))
  if(type=="shared")
  {
    x2$beta <- x1$beta
  } else if(type == "1") {
    x2$beta <- 0
  }
  
  x2 <- x2 %>%
    generate_gwas_ss(nid=nid2, ldobj=ld2)
  
  return(list(x1=x1, x2=x2))
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


sim <- function(lddir1, lddir2, region, nid1, nid2, nsim=100, pshared, pdistinct, p1, hsq1=0.4, window=250000, winnerscurse=TRUE, removecv=TRUE, sim=1, mc.cores=1)
{
  type <- sample(c("shared", "distinct", "1"), nsim, replace=T, prob=c(pshared, pdistinct, p1))
  ld <- organise_ldobj(lddir1, lddir2, region)
  hsq <- 1:nsim / sum(1:nsim) * hsq1
  l <- pblapply(1:nsim, function(i)
  {
    x <- simss(ld[[1]], ld[[2]], nid1, nid2, type[i], hsq[i], window)
    x1 <- x$x1
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
    if(winnerscurse)
    {
      oraw <- prop_overlap(x1raw$bhat, x2raw$bhat, x1raw$se, x2raw$se, 0.05/nsim)
    } else {
      oraw <- prop_overlap(rnorm(nrow(x1raw), x1raw$beta, x1raw$se), x2raw$bhat, x1raw$se, x2raw$se, 0.05/nsim)
    }
    
    fvalsum <- x1$fval + x2$fval
    maxz <- which.max(fvalsum)[1]
    x1maxz <- x1[maxz, ]
    x2maxz <- x2[maxz, ]
    if(winnerscurse)
    {
      omaxz <- prop_overlap(x1maxz$bhat, x2maxz$bhat, x2maxz$se, x2maxz$se, 0.05/nsim)
    } else {
      omaxz <- prop_overlap(rnorm(nrow(x1maxz), x1maxz$beta, x1maxz$se), x2maxz$bhat, x2maxz$se, x2maxz$se, 0.05/nsim)
    }
    return(list(raw=oraw, maxz=omaxz, th1=x1raw$snp, th2=x1maxz$snp))
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
  o1$variants$truth <- type[sig] == "shared"
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
  o2$variants$truth <- type[sig] == "shared"
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
    auc1=o1$auc,
    auc2=o2$auc,
    sim=sim,
    correct_sig1=sum(o1$variants$correct_sig)/length(o1$variants$correct_sig),
    correct_sig2=sum(o2$variants$correct_sig)/length(o2$variants$correct_sig),
    obs_rep1=obs_rep1,
    obs_rep2=obs_rep2,
    exp_rep1=exp_rep1,
    exp_rep2=exp_rep2,
    binom_p1 = binom1$p.value,
    binom_p2 = binom2$p.value,
    fdr1=subset(perf1, !truth)$pow,
    fdr2=subset(perf2, !truth)$pow,
    pow1=subset(perf1, truth)$pow,
    pow2=subset(perf2, truth)$pow,
    prop_same=prop_same
  )
  
  return(res)

}
