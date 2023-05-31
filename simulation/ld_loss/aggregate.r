library(dplyr)
library(jsonlite)
library(here)
config <- read_json(here("simulation/ld_loss/config.json"))
a <- list.files(config$outdir, pattern="result2.rdata", recursive=TRUE)
length(a)

o <- lapply(a, function(f){
    load(file.path(config$outdir, f))
    return(list(snplist, res))
})

res <- lapply(o, function(x) x[[2]]) %>% bind_rows()
snplist <- lapply(o, function(x) x[[1]]) %>% bind_rows()
save(res, snplist, file=file.path(config$resultdir, "results2.rdata"))
