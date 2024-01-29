library(dplyr)
library(tidyr)
library(data.table)
library(R.utils)
library(here)
library(ieugwasr)
library(parallel)
library(GenomicRanges)

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

read_file <- function(m, minmaf=0.01) {
    stopifnot(nrow(m) == 1)
    stopifnot(file.exists(m$fn))
    a <- data.table::fread(m$fn)
    message("Read ", nrow(a), " rows")
    b <- tibble(
        chr = a[[m$chr_col]],
        pos = as.numeric(a[[m$pos_col]]),
        eaf = as.numeric(a[[m$eaf_col]]),
        beta = as.numeric(a[[m$beta_col]]),
        se = as.numeric(a[[m$se_col]]),
        pval = as.numeric(a[[m$pval_col]]),
        ea = a[[m$ea_col]],
        oa = a[[m$oa_col]]
    ) %>% 
    filter(eaf > minmaf & eaf < (1-minmaf)) %>%
    standardise()
    return(b)
}

pool_tophits <- function(rawdat, tophits, metadata, radius = 250000, pthresh = 5e-8, mc.cores = 10) {
    regions <- GRanges(
        seqnames = tophits$chr,
        ranges = IRanges(start=tophits$pos-radius, end=tophits$pos+radius),
        vid=tophits$vid, 
        pop=tophits$pop,
        trait=tophits$trait
    )
    region_list <- lapply(unique(tophits$trait), \(tr) {
        temp <- reduce(subset(regions, trait == tr))
        temp$trait <- tr
        temp
    })

    region_extract <- lapply(1:length(region_list), \(tr) {
        region <- region_list[[tr]]
        mclapply(1:length(region), \(i) {
            message(i, " of ", length(region))
            a <- lapply(1:nrow(metadata), \(j) {
                subset(rawdat[[j]], chr == as.character(seqnames(region)[i]) & pos <= end(region)[i] & pos >= start(region)[i]) %>%
                    mutate(trait = metadata$trait[j], pop = metadata$pop[j], id = metadata$id[j])
            }) %>% bind_rows()
        })
    })

    pool <- lapply(1:length(region_extract), \(tr) {
        region <- region_extract[[tr]]
        mclapply(1:length(region), \(i) {
            target_trait <- region_list[[tr]]$trait[1]
            a <- region[[i]]
            k <- a %>% group_by(vid) %>% summarise(nstudies=n())
            a <- left_join(a, k, by="vid")
            k <- a %>% filter(trait == target_trait) %>%
                group_by(nstudies) %>% 
                summarise(minp = min(pval)) %>% 
                filter(minp < pthresh)
            a <- subset(a, nstudies %in% k$nstudies)
            k <- subset(a, trait == target_trait) %>% 
                mutate(z = abs(beta)/se) %>%
                {subset(., z==max(z))$vid[1]}
            a <- subset(a, vid == k) %>% mutate(target_trait=target_trait)
            return(a)
        }, mc.cores=10) %>% 
            bind_rows()
    }) %>% bind_rows()

    region_list <- lapply(region_list, as_tibble)

    out <- list(region_list=region_list, region_extract=region_extract, tophit_pool=pool)
    return(out)
}

organise_data <- function(metadata,  plink_bin, ld_ref, pthresh=5e-8, minmaf = 0.01, radius = 250000, mc.cores = 1, rawdat=NULL) {
    # read in data

    if(is.null(rawdat)) {
        rawdat <- mclapply(1:nrow(metadata), \(i) read_file(metadata[i,]))
    }

    # get top hits for each
    tophits <- lapply(1:nrow(metadata), \(i) {
        print(i)
        x <- rawdat[[i]] %>% 
            filter(pval < pthresh) %>%
            mutate(rsid = vid)
        if(nrow(x) > 1) {
            ieugwasr::ld_clump(x, plink_bin=plink_bin, bfile=subset(ld_ref, pop == metadata$pop[i])$bfile) %>%
                select(-c(rsid)) %>%
                mutate(pop=metadata$pop[i], trait=metadata$trait[i])
        } else {
            NULL
        }
    }) %>% bind_rows()

    # Get +-250kb region for every tophit
    # Get union of regions
    # Extract regions from each trait
    # Keep SNPs that have at least one GWAS sig and present in all
    # Clump to get tophits
    out <- pool_tophits(rawdat, tophits, metadata, radius = radius, pthresh = pthresh, mc.cores = mc.cores)
    return(out)
}

fixed_effects_meta_analysis_fast <- function(beta_mat, se_mat) {
    w <- 1 / se_mat^2
    beta <- rowSums(beta_mat * w, na.rm=TRUE) / rowSums(w, na.rm=TRUE)
    se <- sqrt(1 / rowSums(w, na.rm=TRUE))
    pval <- pnorm(abs(beta / se), lower.tail = FALSE)
    return(pval)
}

metadata <- readRDS(here("data", "stroke_ldl", "metadata.rds"))

ld_ref <- tibble(
    pop = unique(metadata$pop),
    bfile = here("data", "stroke_ldl", "ref", paste0(pop, "_vid"))
)

# Read LDL in once
metadata_ldl <- subset(metadata, trait == "LDL")
rawdat <- mclapply(1:nrow(metadata_ldl), \(i) read_file(metadata_ldl[i,]))

out <- subset(metadata, what == "outcome")$trait %>%
    unique() %>% {
    lapply(., \(x) {
        message(x)

        temp <- subset(metadata, trait == x)
        rawdat_this <- mclapply(1:nrow(temp), \(i) read_file(temp[i,]))
        rawdat_this <- c(rawdat, rawdat_this)

        a <- organise_data(subset(metadata, trait %in% c("LDL", x)), "plink", ld_ref, mc.cores=10, rawdat=rawdat_this)
        return(a)
    })}

saveRDS(out, file=here("data", "stroke_ldl", "organised_summary_data.rds"))
