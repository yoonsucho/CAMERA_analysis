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
  sim=1:3,
  bxy1=0.2,
  bxy2=0.2,
  mc.cores=16
) %>% filter(
    pshared + pdistinct + p1 == 1
  ) %>%
  mutate(simid=1:n())

dim(param)
o <- mclapply(1:nrow(param), function(i)
  {
    message(i, " of ", nrow(param))
    tryCatch(do.call(sim, args=param[i,]), error=function(e) {print(e); return(NULL)})
  }, mc.cores=1)

oinst <- lapply(o, function(x) x$instruments) %>% bind_rows()
omr <- lapply(o, function(x) x$mr) %>% bind_rows() %>% inner_join(param, by="simid")

save(oinst, omr, file=here("results", "instrument_specificity.rdata"))
