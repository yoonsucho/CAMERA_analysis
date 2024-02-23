library(dplyr)
library(here)

out <- readRDS(here("data", "stroke_ldl", "organised_summary_data.rds"))
names(out)
metadata <- readRDS(here("data", "stroke_ldl", "metadata.rds"))

o <- out[[1]]
inst <- unique(subset(o$tophit_pool, target_trait == "LDL")$vid)
inst_o <- unique(subset(o$tophit_pool, target_trait == "Stroke")$vid)

names(o$region_extract[[1]]) <- inst
names(o$region_extract[[2]]) <- inst_o

instrument_raw <- o$tophit_pool %>% filter(target_trait == "LDL" & trait == "LDL" & nstudies == 5) %>% rename(position="pos", nea="oa", p="pval", rsid="vid")
instrument_outcome <- subset(o$tophit_pool, trait == "Stroke" & target_trait == "LDL" & vid %in% instrument_raw$rsid) %>% rename(position="pos", nea="oa", p="pval", rsid="vid")

# restrict regions to common snps

instrument_raw
instrument_regions <- lapply(unique(instrument_raw$rsid), \(x) {
    a <- o$region_extract[[1]][[x]] %>% 
        filter(trait == "LDL") %>% 
        rename(position="pos", nea="oa", p="pval", rsid="vid") %>%
        group_by(pop) %>% 
        group_split() %>% as.list()
    names(a) <- sapply(a, \(z) z$id[1])
    a
})

instrument_outcome_regions <- lapply(unique(instrument_raw$rsid), \(x) {
    a <- o$region_extract[[1]][[x]] %>% 
        filter(trait == "Stroke") %>% 
        rename(position="pos", nea="oa", p="pval", rsid="vid") %>%
        group_by(pop) %>% 
        group_split() %>% as.list()
    names(a) <- sapply(a, \(z) z$id[1])
    a
})

names(instrument_regions) <- unique(instrument_raw$rsid)
names(instrument_outcome_regions) <- unique(instrument_raw$rsid)

save(instrument_regions, instrument_raw, instrument_outcome, instrument_outcome_regions, file="~/repo/CAMERA/inst/extdata/example-local.rdata")

metadata
metadata2 <- metadata %>% filter(trait %in% c("LDL", "Stroke")) %>% select(what, pop, trait, id, n, ends_with("col"), fn, units) %>%
mutate(fn = basename(fn))
saveRDS(metadata2, file="~/repo/CAMERA/inst/extdata/example-metadata.rds")
