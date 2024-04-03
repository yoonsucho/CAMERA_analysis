library(simulateGP)
library(dplyr)
library(here)
library(jsonlite)

main <- function()
{
  print("HELLOOO")
  ## Input
  config <- jsonlite::read_json("config.json")
  ld_data_dir <- config$ld_data_dir
  region_code <- commandArgs(T)[1]
  outdir <- file.path(config$outdir, region_code)
  dir.create(outdir, recursive=TRUE)
  pops <- c("EUR", "EAS", "SAS", "AFR", "AMR")
  ldobjs <- generate_ldobjs(region_code, outdir, ld_data_dir, pops)
  tempmap <- ldobjs[["EUR"]]$map
  tempmap$ldscore <- colSums(ldobjs$EUR$ld^2)

  load(here("simulation/ld_loss/data","gwashits.rdata"))
  gwashits <- subset(gwashits, code %in% region_code & rsid %in% tempmap$snp)
  matched <- match_maf_ldscore(tempmap, gwashits$rsid) %>%
    left_join(., gwashits, by=c("target"="rsid")) %>%
    rename(rsid=matched)
  snplist <- bind_rows(gwashits, matched)
  print(str(snplist))
  res <- lapply(1:nrow(snplist), function(i) {
    message(i)
    print(snplist[i,] %>% str())
    tophit_cross_ancestry_ld_loss_empirical(tempmap, snplist$rsid[i], snplist$p[i], snplist$n[i], ldobjs)[[2]]
  }) %>% bind_rows() %>% mutate(region=region_code)
  save(snplist, res, file=file.path(outdir, "result.rdata"))
  unlink(file.path(outdir, pops), recursive=TRUE)
}

# For a list of snps find SNPs that match on MAF and LD Score
# Return a data frame of each target SNP and a list of nmatch matching SNPs
match_maf_ldscore <- function(tempmap, snp, bins=10, nmatch=5)
{
  tempmap$ldscore_bin <- cut(tempmap$ldscore, bins)
  tempmap$maf <- tempmap$af
  tempmap$maf[tempmap$maf > 0.5] <- 1 - tempmap$maf[tempmap$maf > 0.5]
  tempmap$maf_bin <- cut(tempmap$maf, bins)
  tempmap$bin <- paste(tempmap$maf_bin, tempmap$ldscore_bin)
  lapply(snp, function(s)
  {
    matched <- tempmap[tempmap$bin %in% tempmap$bin[tempmap$snp==s],]
    tibble(
      target=s,
      matched=matched$snp[sample(1:nrow(matched), min(nrow(matched), nmatch), replace=FALSE)]
    )
  }) %>% bind_rows()
}

generate_ldobjs <- function(region_code, outdir, ld_data_dir, pops = c("EUR", "EAS", "SAS", "AFR", "AMR"))
{
  region <- system.file("extdata/ldetect/EUR.bed", package = "simulateGP") %>% 
    data.table::fread(., header = TRUE) %>% 
    dplyr::mutate(
      chr = as.numeric(gsub("chr", "", chr)), 
      start = as.numeric(start), 
      stop = as.numeric(stop),
      code = paste(chr, start, stop, sep="_")) %>% 
    dplyr::as_tibble() %>%
    dplyr::filter(code == region_code)

  maps <- lapply(pops, function(pop)
  {
    if(!file.exists(file.path(outdir, pop, paste0("ldobj_", region$code, ".rds"))))
    {
      generate_ldobj(file.path(outdir, pop), file.path(ld_data_dir, pop), region)
    } else {
      readRDS(file.path(outdir, pop, "map.rds"))
    }
  })
  names(maps) <- pops
  save(maps, file=file.path(outdir, "maps.rdata"))
  ldobjs <- organise_ldobj(as.list(file.path(outdir, pops, paste0("ldobj_", region$code, ".rds"))), region$code)
  names(ldobjs) <- pops
  glimpse(ldobjs)
  return(ldobjs)
}

subset_ldobj <- function(ldobj, snps)
{
  i <- which(ldobj$map$snp %in% snps)
  ldobj$map <- ldobj$map[i,]
  ldobj$ld <- ldobj$ld[i,i]
  return(ldobj)
}

organise_ldobj <- function(dirlist, region)
{
  ld <- lapply(dirlist, readRDS)
  snplist <- Reduce(intersect, lapply(ld, function(x) x$map$snp))
  ld <- lapply(ld, function(x) subset_ldobj(x, snplist))
  stopifnot(all(ld[[1]]$map$snp == ld[[2]]$map$snp))
  return(ld)
}

b_from_pafn <- function(pval, af, n)
{
  tval <- qnorm(pval, low=F)
  rsq <- tval^2 / (tval^2 + n)
  beta <- sqrt(rsq / (2*af*(1-af)))
  beta
}

prob_y_gt_x <- function(m1, se1, m2, se2)
{
  m <- m1 - m2
  se <- se1^2 + se2^2
  pnorm(0, m, sqrt(se), low=F)
}

tophit_cross_ancestry_ld_loss <- function(tempmap, snp, pval, nEUR, ldobjs)
{
  causalsnp <- which(tempmap$snp == snp)
  tempmap$causal <- 1:nrow(tempmap) == causalsnp
  tempmap$ld <- ldobjs$EUR$ld[causalsnp, ]
  beta <- b_from_pafn(pval, tempmap$af[tempmap$causal], nEUR)
  tempmap$beta <- beta * sign(tempmap$ld) * tempmap$ld^2
  tempmap$se <- expected_se(tempmap$beta, tempmap$af, nEUR, 1)
  # The probability that each SNP has a larger assoc than the causal variant
  tempmap$prob_tophit <- sapply(1:nrow(tempmap), function(i) prob_y_gt_x(tempmap$beta[i], tempmap$se[i], tempmap$beta[causalsnp], tempmap$se[causalsnp]))
  # The scaled probability of being 
  tempmap$prob_tophit2 <- tempmap$prob_tophit / sum(tempmap$prob_tophit)
  r2 <- sapply(ldobjs, function(x){
    x <- subset_ldobj(x, tempmap$snp)$ld[tempmap$causal,]
    sum(x * tempmap$prob_tophit2) / sum(tempmap$prob_tophit2)
  })
  af <- sapply(ldobjs, function(x){
    x <- subset_ldobj(x, tempmap$snp)$map$af[tempmap$causal]
    sum(x * tempmap$prob_tophit2) / sum(tempmap$prob_tophit2)
  })
  return(list(tempmap, tibble(snp=tempmap$snp[tempmap$causal], pop=names(r2), r2=r2, af=af)))
}

tophit_cross_ancestry_ld_loss_empirical_error <- function(tempmap, snp, pval, nEUR, ldobjs, ntry = 10000, seed=sample(1:100000, 1))
{
  set.seed(seed)
  causalsnp <- which(tempmap$snp == snp)
  tempmap$causal <- 1:nrow(tempmap) == causalsnp
  tempmap$ld <- ldobjs$EUR$ld[causalsnp, ]
  beta <- b_from_pafn(pval, tempmap$af[tempmap$causal], nEUR)

  tempmap$beta2 <- beta * sign(tempmap$ld) * tempmap$ld^2
  tempmap$se <- expected_se(tempmap$beta, tempmap$af, nEUR, 1)
  out <- sapply(1:ntry, function(i)
  {
    tempmap$bhat <- sample_beta(tempmap$beta, tempmap$se)
    tempmap$fval <- (tempmap$bhat / tempmap$se)^2
    tempmap$snp[which.max(tempmap$fval)[1]]
  })
  tempmap <- prop.table(table(out)) %>% as_tibble %>% left_join(tempmap, ., by=c("snp"="out"))
  tempmap$n[is.na(tempmap$n)] <- 0
  tempmap$prob_tophit2 <- tempmap$n
  r2 <- sapply(ldobjs, function(x){
    x <- subset_ldobj(x, tempmap$snp)$ld[tempmap$causal,]
    sum(x * tempmap$prob_tophit2) / sum(tempmap$prob_tophit2)
  })
  af <- sapply(ldobjs, function(x){
    x <- subset_ldobj(x, tempmap$snp)$map$af[tempmap$causal]
    sum(x * tempmap$prob_tophit2) / sum(tempmap$prob_tophit2)
  })
  frac_causal <- sapply(ldobjs, function(x){
    x <- subset_ldobj(x, tempmap$snp)$prob[tempmap$causal,]
  })
  # r2 = average rsq of the detected hits with the causal variant
  # af = average allele frequency of the detected hits
  # n_nonzero = number of variants with non-zero chance of being selected
  # causal_prob = probability of the causal variant being selected
  tib <- tibble(snp=tempmap$snp[tempmap$causal], pop=names(r2), pval=pval, n=nEUR, seed=seed, r2=r2, af=af, n_nonzero = NA, causal_prob = NA)
  tib$n_nonzero[tib$pop == "EUR"] <- sum(tempmap$prob_tophit2 != 0)
  tib$causal_prob[tib$pop == "EUR"] <- tempmap$prob_tophit2[tempmap$causal]
  return(list(tempmap, tib))
}

tophit_cross_ancestry_ld_loss_empirical <- function(tempmap, snp, pval, nEUR, ldobjs, ntry = 10000, seed=sample(1:100000, 1), max_snps=2000)
{
  print(str(tempmap))
  set.seed(seed)
  causalsnp <- which(tempmap$snp == snp)
  print(causalsnp)
  tempmap$causal <- 1:nrow(tempmap) == causalsnp
  min_coord <- max(1, round(causalsnp - (max_snps/2)))
  max_coord <- min(nrow(tempmap), round(causalsnp + (max_snps/2)))
  message(min_coord, "-", max_coord)
  tempmap <- tempmap[min_coord:max_coord,]
  causalsnp <- which(tempmap$causal)
  message("Number of SNPs: ", nrow(tempmap))
  # tempmap$ld <- ldobjs$EUR$ld[causalsnp, ]
  beta <- b_from_pafn(pval, tempmap$af[tempmap$causal], nEUR)
  tempmap$beta <- 0
  tempmap$beta[tempmap$causal] <- beta
  xvar <- tempmap$af * (1-tempmap$af) * 2
  eurld <- subset_ldobj(ldobjs$EUR, tempmap$snp)
  tempmap$beta_ld <- (diag(1/xvar) %*% eurld$ld %*% diag(xvar) %*% tempmap$beta) %>% drop()
  tempmap$ld <- eurld$ld[causalsnp, ]
  #tempmap$beta_ld <- beta * sign(tempmap$ld) * tempmap$ld^2
  tempmap$se <- expected_se(tempmap$beta_ld, tempmap$af, nEUR, 1)
  semat <- diag(tempmap$se) %*% eurld$ld %*% diag(tempmap$se)
  bhat <- MASS::mvrnorm(ntry, mu=tempmap$beta_ld, Sigma=semat)
  fval <- (t(bhat) / tempmap$se)^2
  topsnp <- tempmap$snp[apply(fval, 2, function(x) which.max(x))]
  tempmap$ord <- 1:nrow(tempmap)
  tempmap <- prop.table(table(topsnp)) %>% as_tibble %>% dplyr::select(snp=topsnp, prob_tophit=n) %>% left_join(tempmap, .) %>% arrange(ord)
  tempmap$prob_tophit[is.na(tempmap$prob_tophit)] <- 0
  r2 <- sapply(ldobjs, function(x){
    x <- subset_ldobj(x, tempmap$snp)$ld[tempmap$causal,]^2
    sum(x * tempmap$prob_tophit) / sum(tempmap$prob_tophit)
  })
  af <- sapply(ldobjs, function(x){
    x <- subset_ldobj(x, tempmap$snp)$map$af[tempmap$causal]
    sum(x * tempmap$prob_tophit) / sum(tempmap$prob_tophit)
  })
  frac_causal <- sapply(ldobjs, function(x){
    x <- subset_ldobj(x, tempmap$snp)$prob[tempmap$causal,]
  })
  # r2 = average rsq of the detected hits with the causal variant
  # af = average allele frequency of the detected hits
  # n_nonzero = number of variants with non-zero chance of being selected
  # causal_prob = probability of the causal variant being selected
  tib <- tibble(snp=tempmap$snp[tempmap$causal], pop=names(r2), pval=pval, n=nEUR, seed=seed, r2=r2, af=af, n_nonzero = NA, causal_prob = NA)
  tib$n_nonzero[tib$pop == "EUR"] <- sum(tempmap$prob_tophit != 0)
  tib$causal_prob[tib$pop == "EUR"] <- tempmap$prob_tophit[tempmap$causal]
  return(list(tempmap, tib))
}

main()
