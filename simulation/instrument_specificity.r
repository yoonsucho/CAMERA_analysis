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
  nid2=c(50000, 100000, 200000),
  pshared=c(0, 0.5),
  pdistinct=c(0, 0.5),
  p1=0.5,
  nsim=100,
  hsq1=0.4,
  window=250000,
  sim=1:5
) %>% filter(
  pshared + pdistinct + p1 == 1
)

o <- mclapply(1:nrow(param), function(i)
  {
    message(i, " of ", nrow(param))
    do.call(sim, args=param[i,])
  }, mc.cores=16) %>% bind_rows()

save(o, file=here("results", "instrument_specificity.rdata"))
