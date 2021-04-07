library(tidyverse)
library(ieugwasr)
library(simulateGP)
library(lavaan)
library(TwoSampleMR)
library(semPlot)
library(tidySEM)
library(MASS)
library(dplyr)
library(parallel)

##------------------------------original function------------------------------##
sim1 <- function(nsnp, hsq1, hsq2, rg, nid1, nid2, maf1, maf2, biv1, biv2, inst_pval1, inst_pval2)
{
  # Generate effects for exposures in pop1 and pop2
  s1 <- sqrt(hsq1 / nsnp)
  s2 <- sqrt(hsq2 / nsnp)
  sigma <- matrix(c(1, rg, rg, 1), 2, 2)
  eff <- mvrnorm(nsnp, c(0,0), sigma)
  
  gx1 <- eff[,1] * s1
  gx2 <- eff[,2] * s2

  # create genotypes
  g1 <- make_geno(nid1, nsnp, maf1)
  g2 <- make_geno(nid2, nsnp, maf2)
  
  # create confounder
  u1 <- rnorm(nid1)
  u2 <- rnorm(nid2)
  
  # create x
  x1 <- make_phen(c(0.1, gx1), cbind(u1, g1))
  x2 <- make_phen(c(0.1, gx2), cbind(u2, g2))

  # create y
  y1 <- make_phen(c(0.1, biv1), cbind(u1, x1))
  y2 <- make_phen(c(0.1, biv2), cbind(u2, x2))

  # Summary data for pop1 
  ss1 <- get_effs(x1, y1, g1)
  
  # Summary data for pop2
  ss2 <- get_effs(x2, y2, g2)

  # Organise the data
  dat <- tibble(
  	y1 = ss1$beta.outcome,
  	y2 = ss2$beta.outcome,
  	x1 = ss1$beta.exposure,
  	x2 = ss2$beta.exposure
  )
  
  index <- ss1$pval.exposure < inst_pval1 & ss2$pval.exposure < inst_pval2
  dat <- dat[index,]

  model1 <- '
  y1 ~ biv*x1
  y2 ~ biv*x2
  '
  mod1 <- sem(model1, data=dat)

  model2 <- '
  y1 ~ biv_1*x1
  y2 ~ biv_2*x2
  '
  mod2 <- sem(model2, data=dat)
  
  out <- list()
  invisible(capture.output(s1 <- summary(mod1, fit.measures=TRUE)))
  invisible(capture.output(s2 <- summary(mod2, fit.measures=TRUE)))
  mod1 <- s1$FIT['aic']
  mod2 <- s2$FIT['aic']
  out$aic1 <- s1$FIT['aic']
  out$aic2 <- s2$FIT['aic']
  out$aic_diff <- mod1 - mod2

  out$mod1_eff1 <- s1$PE$est[1]
  out$mod1_eff2 <- s1$PE$est[2]
  out$mod2_eff1 <- s2$PE$est[1]
  out$mod2_eff2 <- s2$PE$est[2]

  out$mod1_pval1 <- s1$PE$pvalue[1]
  out$mod1_pval2 <- s1$PE$pvalue[2]
  out$mod2_pval1 <- s2$PE$pvalue[1]
  out$mod2_pval2 <- s2$PE$pvalue[2]

  return(as_tibble(out))
}
##----------------------------------------------------------------------------##




sim_params <- expand.grid(
  nid1 = c(10000, 20000),
  nid2 = c(10000, 20000),
  nsnp = 100,
  hsq1 = c(0.1, 0.4),
  hsq2 = c(0.1, 0.4),
  maf1 = c(0.1, 0.5),
  maf2 = c(0.1, 0.5),
  biv1 = seq(-0.1, 0.1, by=0.01),
  biv2 = seq(-0.1, 0.1, by=0.01),
  inst_pval1 = c(1, 0.05, 0.05/100),
  inst_pval2 = c(1, 0.05, 0.05/100),
  rg = c(0, 0.5, 1),
  sim = c(1:1000)
)

sim_params$sim <- 1:nrow(sim_params)
args <- commandArgs(T)
chunk <- as.numeric(args[1])
chunksize <- as.numeric(args[2])
ncores <- as.numeric(args[3])
a <- (chunk - 1) * chunksize + 1
b <- min(chunk * chunksize, nrow(sim_params))
message("Running ", a, " to ", b)
sim_params <- sim_params[a:b, ]



#dim(sim_params)

# l <- list()
# for(i in 1:nrow(sim_params))
# {
#   res <- 
#   l[[i]] <- bind_cols(sim_params[i,], res)
# }
# l <- bind_rows(l)


res <-  mclapply(1:nrow(sim_params), function(i)
                {
                  set.seed(sim_params$sim[i])
                  bind_cols(sim_params[i,], do.call(sim1, args=as.list(sim_params[i, ])))
                }, mc.cores=ncores) %>% bind_rows()
                

# 1. Are the effect estimates different? - (AIC difference)
  # - how different do the effects need to be for AIC < -2
# 2. Are the effect estimates unbiased? - (estimate - true)
  # - Under what scenarios are the estimates wrong
# 3. Are the effect estimates significant - (pval < 0.05)
  # - Under what scenarios are FDR high or power low

save(res, param, file=paste0("/work/yc16575/mr.trans/result/sim_", chunk, ".rdata"))
