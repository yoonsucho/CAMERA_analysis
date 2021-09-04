library(tidyverse)
library(here)
nchunk <- as.numeric(commandArgs(T))

l1 <- list()
l2 <- list()
for(i in 1:nchunk)
{
	load(here("results", "instrument_specificity", paste0(i, ".rdata")))
	l1[[i]] <- oinst
	l2[[i]] <- omr
}

oinst <- bind_rows(l1)
omr <- bind_rows(l2)
save(oinst, omr, file=here("results", "instrument_specificity.rdata"))
