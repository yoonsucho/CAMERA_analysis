match_maf_ldscore <- function(tempmap, snp, bins=10, nmatch=5) {
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

generate_ldobjs <- function(region_code, outdir, ld_data_dir, pops = c("EUR", "EAS", "SAS", "AFR", "AMR")) {
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

subset_ldobj <- function(ldobj, snps) {
  i <- which(ldobj$map$snp %in% snps)
  ldobj$map <- ldobj$map[i,]
  ldobj$ld <- ldobj$ld[i,i]
  return(ldobj)
}

organise_ldobj <- function(dirlist, region) {
  ld <- lapply(dirlist, readRDS)
  snplist <- Reduce(intersect, lapply(ld, function(x) x$map$snp))
  ld <- lapply(ld, function(x) subset_ldobj(x, snplist))
  stopifnot(all(ld[[1]]$map$snp == ld[[2]]$map$snp))
  return(ld)
}

b_from_pafn <- function(pval, af, n) {
  tval <- qnorm(pval, low=F)
  rsq <- tval^2 / (tval^2 + n)
  beta <- sqrt(rsq / (2*af*(1-af)))
  beta
}

prob_y_gt_x <- function(m1, se1, m2, se2) {
  m <- m1 - m2
  se <- se1^2 + se2^2
  pnorm(0, m, sqrt(se), low=F)
}

tophit_cross_ancestry_ld_loss <- function(tempmap, snp, pval, nEUR, ldobjs) {
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

tophit_cross_ancestry_ld_loss_empirical_error <- function(tempmap, snp, pval, nEUR, ldobjs, ntry = 10000, seed=sample(1:100000, 1)) {
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

tophit_cross_ancestry_ld_loss_empirical <- function(tempmap, snp, pval, nEUR, ldobjs, ntry = 10000, seed=sample(1:100000, 1), max_snps=2000) {
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

#' Fixed effects meta analysis vectorised across multiple SNPs
#' 
#' @param beta_mat Matrix of betas - rows are SNPs, columns are studies
#' @param se_mat Matrix of SEs - rows are SNPs, columns are studies
#' 
#' @return list of meta analysis betas and SEs
fema <- function(beta_mat, se_mat) {
  w <- 1/se_mat^2
  beta <- rowSums(beta_mat * w) / rowSums(w)
  se <- sqrt(1/rowSums(w))
  return(list(beta=beta, se=se))
}
# betas_mat <- rbind(c(0.4, 0.6, 0.8, 1.0, 1.2), sample(c(0.4, 0.6, 0.8, 1.0, 1.2)))
# se_mat <- rbind(c(0.2, 0.1, 0.3, 0.4, 0.2), sample(c(0.2, 0.1, 0.3, 0.4, 0.2)))

# fema(betas_mat, se_mat)
# rma(yi=betas_mat[1,], sei=se_mat[1,], method="FE")
# rma(yi=betas_mat[2,], sei=se_mat[2,], method="FE")

# fixed_effects_meta(betas, se)


generate_tophit_probs <- function(tempmap, snp, pval, nEUR, ldobjs, ntry = 10000, seed=sample(1:100000, 1), max_snps=2000){
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
  semat <- diag(tempmap$se) %*% eurld$ld %*% diag(tempmap$se)
  bhat <- MASS::mvrnorm(ntry, mu=tempmap$beta_ld, Sigma=semat)
  fval <- (t(bhat) / tempmap$se)^2
  topsnp <- tempmap$snp[apply(fval, 2, function(x) which.max(x))]
  tempmap$ord <- 1:nrow(tempmap)
  tempmap <- prop.table(table(topsnp)) %>% as_tibble %>% dplyr::select(snp=topsnp, prob_tophit=n) %>% left_join(tempmap, .) %>% arrange(ord)
  tempmap$prob_tophit[is.na(tempmap$prob_tophit)] <- 0
  return(tempmap)
}


generate_tophit_probs_ma <- function(tempmap, snp, pval, nEUR, ldobjs, ntry = 10000, seed=sample(1:100000, 1), max_snps=2000, prop_afr){
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
  tempmap$se <- expected_se(tempmap$beta_ld, tempmap$af, nEUR * (1-prop_afr), 1)

  # add x% african sample 
  afrld <- subset_ldobj(ldobjs$AFR, tempmap$snp)
  xvar_afr <- afrld$map$af * (1-afrld$map$af) * 2
  tempmap$beta_ld_afr <- (diag(1/xvar_afr) %*% afrld$ld %*% diag(xvar_afr) %*% tempmap$beta) %>% drop()
  tempmap$ld_afr <- afrld$ld[causalsnp, ]
  #tempmap$beta_ld <- beta * sign(tempmap$ld) * tempmap$ld^2
  tempmap$se_afr <- expected_se(tempmap$beta_ld_afr, afrld$map$af, nEUR * prop_afr, 1)
  semat_eur <- diag(tempmap$se) %*% eurld$ld %*% diag(tempmap$se)
  semat_afr <- diag(tempmap$se_afr) %*% afrld$ld %*% diag(tempmap$se_afr)
  bhat_eur <- MASS::mvrnorm(ntry, mu=tempmap$beta_ld, Sigma=semat_eur)
  bhat_afr <- MASS::mvrnorm(ntry, mu=tempmap$beta_ld_afr, Sigma=semat_afr)

  fval <- sapply(1:ntry, function(i){
    o <- fema(cbind(bhat_eur[i, ], bhat_afr[i,]), cbind(tempmap$se, tempmap$se_afr))
    (o$beta/o$se)^2
  })
  topsnp <- tempmap$snp[apply(fval, 2, function(x) which.max(x))]
  tempmap$ord <- 1:nrow(tempmap)
  tempmap <- prop.table(table(topsnp)) %>% as_tibble %>% dplyr::select(snp=topsnp, prob_tophit=n) %>% left_join(tempmap, .) %>% arrange(ord)
  tempmap$prob_tophit[is.na(tempmap$prob_tophit)] <- 0
  return(tempmap)
}

eval_tophit_probs <- function(ldobjs, tm) {
  r2 <- sapply(ldobjs, function(x){
    x <- subset_ldobj(x, tm$snp)$ld[tm$causal,]^2
    sum(x * tm$prob_tophit) / sum(tm$prob_tophit)
  })
  af <- sapply(ldobjs, function(x){
    x <- subset_ldobj(x, tm$snp)$map$af[tm$causal]
    sum(x * tm$prob_tophit) / sum(tm$prob_tophit)
  })
  frac_causal <- sapply(ldobjs, function(x){
    x <- subset_ldobj(x, tm$snp)$prob[tm$causal,]
  })
  # r2 = average rsq of the detected hits with the causal variant
  # af = average allele frequency of the detected hits
  # n_nonzero = number of variants with non-zero chance of being selected
  # causal_prob = probability of the causal variant being selected
  tib <- tibble(snp=tm$snp[tm$causal], pop=names(r2), r2=r2, af=af, n_nonzero = NA, causal_prob = NA)
  tib$n_nonzero[tib$pop == "EUR"] <- sum(tm$prob_tophit != 0)
  tib$causal_prob[tib$pop == "EUR"] <- tm$prob_tophit[tm$causal]
  return(tib)
}
