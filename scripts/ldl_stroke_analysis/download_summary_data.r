library(dplyr)
library(yaml)
library(glue)
library(here)

# Get list of stroke traits
bn <- glue("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST901045{35:63}/GCST901045{35:63}_buildGRCh37.tsv.gz-meta.yaml")

d <- lapply(bn, read_yaml)

t <- lapply(d, \(x) {
    tibble(
        what = "outcome",
        id_gc = x$gwas_id,
        trait = x$trait_description,
        pop = x$samples[[1]]$sample_ancestry_category,
        id = paste(trait, pop),
        n = x$samples[[1]]$sample_size,
        nsample = length(x$samples),
        chr_col = 1,
        pos_col = 2,
        eaf_col = 3,
        beta_col = 4,
        se_col = 5,
        pval_col = 6, 
        ea_col = 9,
        oa_col = 10,
        units = "logOR"
    )
}) %>% bind_rows() %>%
mutate(
    fn_yaml = bn, 
    url=gsub("-meta.yaml", "", fn_yaml), 
    fn=here("data", "stroke_ldl", "raw", basename(url))) %>%
filter(nsample == 1, trait == "Stroke")

# Manual step
t$pop <- c("EUR", "EAS", "AFR", "AMR", "SAS")

dir.create(here("data", "stroke_ldl", "raw"), recursive=T)

lapply(1:nrow(t), \(i) {
    message(t$fn[i])
    download.file(t$url[i], destfile = t$fn[i])
})


# LDL traits
ldl <- tibble(
    what = "exposure",
    url = c(
        "https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/LDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",
        "https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/LDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz",
        "https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",
        "https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/LDL_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz",
        "https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/LDL_INV_SAS_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
    pop = c("AFR", "EAS", "EUR", "AMR", "SAS"),
    trait = "LDL",
    id = paste(trait, pop),
    n = c(91144, 79831, 900191, 44790, 30412),
    rsid_col = 1,
    chr_col = 2,
    pos_col = 3,
    eaf_col = 8,
    beta_col = 9,
    se_col = 10,
    pval_col = 12, 
    ea_col = 5,
    oa_col = 4,
    fn = here("data", "stroke_ldl", "raw", basename(url)),
    units = "SD"
)

options(timeout=400)
lapply(ldl$url[1], \(x) download.file(url=x, destfile = here("data", "stroke_ldl", "raw", basename(x))))


# Get sample sizes for LDL - a bit post hoc
# for i in ../../data/stroke_ldl/raw/LDL_INV_*
# do
#     zcat $i | head -n 1000 | awk '{ print $6 }' | sort -n | tail -n 1
# done


metadata <- bind_rows(ldl, t)
metadata
saveRDS(metadata, file=here("data", "stroke_ldl", "metadata.rds"))
