library(here)
library(tidyverse)
source(here("simulation", "instrument_specificity_functions.r"))
load(here("data", "ld", "maps.rdata"))

regions <- unique(mapEUR$region)

param <- expand.grid(
  lddir1=here("data", "ld", "ldEUR"),
  lddir2=here("data", "ld", c("ldEUR", "ldEAS", "ldAFR")),
  region=regions,
  nid1=100000,
  nid2=c(12500, 25000, 50000, 100000, 200000),
  removecv=c(TRUE, FALSE),
  pshared=c(0, 0.5),
  pdistinct=c(0, 0.5),
  p1=0.5,
  nsim=100,
  hsq1=0.1,
  window=250000,
  sim=1:10,
  bxy1=seq(0, 0.05, by=0.01),
  bxy2=seq(0, 0.05, by=0.01),
  mc.cores=16
) %>% filter(
    pshared + pdistinct + p1 == 1
  ) %>%
  mutate(simid=1:n())


args <- commandArgs(T)
chunk <- as.numeric(args[1]) # starts from 0
chunksize <- as.numeric(args[2])
nchunk <- ceiling(nrow(param) / chunksize)
start <- chunk * chunksize + 1
end <- min((chunk + 1) * chunksize, nrow(param))

message("total:", nrow(param))
message("start: ", start)
message("end: ", end)

param <- param[start:end,]

output <- here("results", "instrument_specificity", paste0(chunk, ".rdata"))
output_int <- paste0(output, ".int")

if(file.exists(output))
{
  message("already complete")
  q()
}

if(file.exists(output_int))
{
  message("Previous run already exists")
  load(output_int)
  a <- max(sapply(oinst1, function(x) x$simid))
  message(sum(param$simid > a), " out of ", nrow(param), " remaining")
  j <- which(param$simid > a)[1]
} else {
  message("New run")
  oinst1 <- list()
  omr1 <- list()
  j <- 1
}

for(i in j:nrow(param))
{
  message(i, " of ", nrow(param))
  res <- tryCatch(do.call(sim, args=param[i,]), error=function(e) {print(e); return(NULL)})
  oinst1[[i]] <- res$instruments
  oinst1[[i]]$simid <- param$simid[i]
  omr1[[i]] <- res$mr
  omr1[[i]]$simid <- param$simid[i]
  save(oinst1, omr1, file=output_int)
}



oinst <- bind_rows(oinst1)
omr <- bind_rows(omr1) %>% inner_join(param, by="simid")

save(oinst, omr, file=output)

