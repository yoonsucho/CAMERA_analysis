library(simulateGP)

set.seed(12345)
n_regions <- 10
regions <- system.file(paste0("extdata/ldetect/EUR.bed"), package = "simulateGP") %>% 
        data.table::fread(., header = TRUE) %>% dplyr::mutate(chr = as.numeric(gsub("chr", 
        "", chr)), start = as.numeric(start), stop = as.numeric(stop)) %>% 
        dplyr::as_tibble()
regions <- regions[sample(1:nrow(regions), 10, replace=FALSE), ]

bEUR <- "../data/ld/EUR"
bAFR <- "../data/ld/AFR"
bEAS <- "../data/ld/EAS"

varrefEUR <- simulateGP:::variant_reference(bEUR)
varrefAFR <- simulateGP:::variant_reference(bAFR)
varrefEAS <- simulateGP:::variant_reference(bEAS)

mapEUR <- generate_ldobj("../data/ld/ldEUR", bEUR, regions)
mapEAS <- generate_ldobj("../data/ld/ldEAS", bEAS, regions)
mapAFR <- generate_ldobj("../data/ld/ldAFR", bAFR, regions)

save(mapEUR, mapEAS, mapAFR, file="../data/ld/maps.rdata")
