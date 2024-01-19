library(dplyr)
library(data.table)
library(R.utils)
library(here)


# Instrument extraction


# E_i = exposure GWAS for population i
# O_i = exposure GWAS for population i

# T(E_i)_j = Tophits variant j for GWAS E_i for population i

# T(E_i)_ij = Association in exposure GWAS for population i, for tophit variant j for population i

# Outcome tophit pool in all outcomes

# Exposure tophit pool in all exposures
# Exposure tophit pool in all outcomes
# Exposure region pool in all regions




standardise <- function(d, ea_col="ea", oa_col="oa", beta_col="beta", eaf_col="eaf", chr_col="chr", pos_col="pos", vid_col="vid") {
    toflip <- d[[ea_col]] > d[[oa_col]]
    d[[eaf_col]][toflip] <- 1 - d[[eaf_col]][toflip]
    d[[beta_col]][toflip] <- d[[beta_col]][toflip] * -1
    temp <- d[[oa_col]][toflip]
    d[[oa_col]][toflip] <- d[[ea_col]][toflip]
    d[[ea_col]][toflip] <- temp
    d[[vid_col]] <- paste0(d[[chr_col]], ":", d[[pos_col]], "_", d[[ea_col]], "_", d[[oa_col]])
    d
}

read_file <- function(metadata) {
    stopifnot(nrow(metadata) == 1)
    a <- data.table::fread(metadata$fn, nrows=1000)
    message("Read ", nrow(a), " rows")
    b <- tibble(
        chr = a[[metadata$chr_col]],
        pos = a[[metadata$pos_col]],
        eaf = a[[metadata$eaf_col]],
        beta = a[[metadata$beta_col]],
        se = a[[metadata$se_col]],
        pval = a[[metadata$pval_col]],
        ea = a[[metadata$ea_col]],
        oa = a[[metadata$oa_col]]
    ) %>% standardise()
    return(b)
}



generate_dataset <- function(metadata, meta, plink_bin, bfile, pthresh=5e-8) {

    rawdat <- lapply(1:nrow(metadata), \(i) read_file(metadata[i,]))

    for(i in 1:nrow(metadata)) {
        
    }
}

metadata <- readRDS(here("data", "stroke_ldl", "metadata.rds"))


