set.seed(12345)
library(simulateGP)
library(dplyr)
library(here)
library(jsonlite)

main <- function()
{
  ## Input
  config <- jsonlite::read_json("config.json")
  ld_data_dir <- config$ld_data_dir
  region_id <- 1
  outdir <- file.path(config$outdir, region)

  pops <- c("EUR", "EAS", "SAS", "AFR", "AMR")
  pvals <- 10^-seq(8, 64, by=4)

  ldobjs <- generate_ldobjs(region_id, outdir, ld_data_dir, pops)
  tempmap <- ldobjs[["EUR"]]$map
  tempmap$ldscore <- colSums(ldobjs$EUR$ld)

  res <- lapply(pvals, function(pval)
  {
    message(pval)
    lapply(tempmap$snp, function(snp)
    {
      tophit_cross_ancestry_ld_loss(tempmap, snp, pval, 100000, ldobjs)[[2]]
    }) %>% bind_rows() %>% mutate(pval=pval)
  }) %>% bind_rows()
  save(res, file=file.path(outdir, "result.rdata"))
}

generate_ldobjs <- function(region_id, outdir, ld_data_dir, pops = c("EUR", "EAS", "SAS", "AFR", "AMR"))
{
  region <- system.file("extdata/ldetect/EUR.bed", package = "simulateGP") %>% 
    data.table::fread(., header = TRUE) %>% 
    dplyr::mutate(chr = as.numeric(gsub("chr", "", chr)), start = as.numeric(start), stop = as.numeric(stop)) %>% 
    dplyr::as_tibble() %>%
    {.[region_id,]}

  maps <- lapply(pops, function(pop) generate_ldobj(file.path(outdir, pop), file.path(ld_data_dir, pop), region))
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



main()