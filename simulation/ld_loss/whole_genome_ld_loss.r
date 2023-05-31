# here::i_am("simulation/ld_loss/whole_genome_ld_loss.r")
library(here)
library(simulateGP)
library(dplyr)
library(jsonlite)
source(here("simulation/ld_loss/whole_genome_ld_loss_functions.r"))

main <- function()
{
  ## Input
  config <- jsonlite::read_json("config.json")
  ld_data_dir <- config$ld_data_dir
  # region_code <- "12_106958748_109025901"
  region_code <- commandArgs(T)[1]
  message(region_code)
  outdir <- here("simulation", "ld_loss", config$outdir, region_code)
  if(file.exists(file.path(outdir, "result2.rdata"))) 
  {
    message("file already exists")
    q()
  }
  dir.create(outdir, recursive=TRUE)
  pops <- c("EUR", "EAS", "SAS", "AFR", "AMR")
  ldobjs <- generate_ldobjs(region_code, outdir, ld_data_dir, pops)
  tempmap <- ldobjs[["EUR"]]$map
  tempmap$ldscore <- colSums(ldobjs$EUR$ld^2, na.rm=T)

  load(here("simulation", "ld_loss", "data","gwashits.rdata"))
  gwashits <- subset(gwashits, code %in% region_code & rsid %in% tempmap$snp)
  matched <- match_maf_ldscore(tempmap, gwashits$rsid) %>%
    left_join(., gwashits, by=c("target"="rsid")) %>%
    rename(rsid=matched)
  snplist <- bind_rows(gwashits, matched)
  print(str(snplist))
  # i <- 1
  res <- lapply(1:nrow(snplist), function(i) {
    message("RUNNING", i, " of ", nrow(snplist))
    o <- tophit_cross_ancestry_ld_loss_empirical(tempmap, snplist$rsid[i], snplist$p[i], snplist$n[i], ldobjs)
    tempmap1 <- generate_tophit_probs(tempmap, snplist$rsid[i], snplist$p[i], snplist$n[i], ldobjs)
    tempmap2 <- generate_tophit_probs_ma(tempmap, snplist$rsid[i], snplist$p[i], snplist$n[i], ldobjs, prop_afr=0.1)
    tempmap3 <- generate_tophit_probs_ma(tempmap, snplist$rsid[i], snplist$p[i], snplist$n[i], ldobjs, prop_afr=0.5)
    bind_rows(
      eval_tophit_probs(ldobjs, tempmap1) %>% mutate(prop_afr=0, pval=snplist$p[i], n=snplist$n[i]),
      eval_tophit_probs(ldobjs, tempmap2) %>% mutate(prop_afr=0.1, pval=snplist$p[i], n=snplist$n[i]),
      eval_tophit_probs(ldobjs, tempmap3) %>% mutate(prop_afr=0.5, pval=snplist$p[i], n=snplist$n[i])
    )
  }) %>% bind_rows() %>% mutate(region=region_code)
  save(snplist, res, file=file.path(outdir, "result2.rdata"))
  unlink(file.path(outdir, pops), recursive=TRUE)
}

# For a list of snps find SNPs that match on MAF and LD Score
# Return a data frame of each target SNP and a list of nmatch matching SNPs


main()
