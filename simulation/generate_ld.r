set.seed(12345)

library(simulateGP)
library(here)

n_regions <- 10
regions <- system.file("extdata/ldetect/EUR.bed", package = "simulateGP") %>% 
        data.table::fread(., header = TRUE) %>% dplyr::mutate(chr = as.numeric(gsub("chr", 
        "", chr)), start = as.numeric(start), stop = as.numeric(stop)) %>% 
        dplyr::as_tibble()
regions <- regions[sample(1:nrow(regions), 10, replace=FALSE), ]

mapEUR <- generate_ldobj(here("data", "ld", "ldEUR"), here("data", "ld", "EUR"), regions)
mapEAS <- generate_ldobj(here("data", "ld", "ldEAS"), here("data", "ld", "EAS"), regions)
mapAFR <- generate_ldobj(here("data", "ld", "ldAFR"), here("data", "ld", "AFR"), regions)
mapSAS <- generate_ldobj(here("data", "ld", "ldSAS"), here("data", "ld", "SAS"), regions)
mapAMR <- generate_ldobj(here("data", "ld", "ldAMR"), here("data", "ld", "AMR"), regions)

save(mapEUR, mapEAS, mapAFR, mapSAS, mapAMR, file=here("data", "ld", "maps.rdata"))

outdir <- here("data", "ld", "ldEUR")
bfile <- here("data", "ld", "EUR")

generate_ldobj <- function (outdir, bfile, regions, plink_bin = genetics.binaRies::get_plink_binary(), nthreads = 1)
{
    dir.create(outdir)
    codes <- paste0(gsub("chr", "", regions$chr), "_", regions$start,
        "_", regions$stop)
    cmd <- paste0("cat ", bfile, ".fam | wc -l")
    nref <- system(cmd, intern = TRUE) %>% trimws() %>% as.numeric()
    map <- parallel::mclapply(1:nrow(regions), function(i) {
        message(i, " of ", nrow(regions))
        out <- get_ld(chr = gsub("chr", "", regions$chr[i]),
            from = regions$start[i], to = regions$stop[i], bfile = bfile,
            plink_bin = plink_bin, nref = nref)
        if (!is.null(out)) {
            fn <- file.path(outdir, paste0("ldobj_", codes[i],
                ".rds"))
            saveRDS(out, file = fn, compress = TRUE)
            out$map$region <- codes[i]
            return(out$map)
        }
        else {
            return(NULL)
        }
    }, mc.cores = nthreads) %>% dplyr::bind_rows()
    print(dim(map))
    saveRDS(map, file.path(outdir, "map.rds"))
    return(map)
}
