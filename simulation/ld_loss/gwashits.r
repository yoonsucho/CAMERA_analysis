library(ieugwasr)
library(dplyr)
library(simulateGP)
library(here)

# Get all GWAS hits for UKB-b and UKB-d
a <- ieugwasr::gwasinfo()
idlist <- a %>%
    filter(grepl("ukb-b", .$id) | grepl("ukb-d", .$id)) %>%
    {.$id}

gwashits <- tophits(idlist)
summary(gwashits$eaf)

# Need to map these to the LD regions
region <- system.file("extdata/ldetect/EUR.bed", package = "simulateGP") %>% 
    data.table::fread(., header = TRUE) %>% 
    dplyr::mutate(chr = as.numeric(gsub("chr", "", chr)), start = as.numeric(start), stop = as.numeric(stop)) %>% 
    dplyr::as_tibble() %>%
    mutate(code = paste0(chr, "_", start, "_", stop))
str(region)

gwashits <- gwashits %>%
    group_by(chr) %>%
    do({
        x <- .
        reg <- subset(region, chr == x$chr[1])
        for(i in 1:nrow(x))
        {
            j <- which(reg$start <= x$position[i] & reg$stop > x$position[i])
            x$code[i] <- reg$code[j]
        }
        x
    })

save(gwashits, file=here("simulation/ld_loss/data","gwashits.rdata"))
write.table(region, file=here("simulation/ld_loss/data", "region.txt"), row=F, col=F, qu=F)
