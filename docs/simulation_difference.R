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
  
  return(list(model1=mod1, model2=mod2))
}
##----------------------------------------------------------------------------##


#--------------------------sample size--------------------------#
#sample sizes are equal

sample1 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#pop1 has less samples (75%)
sample2 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=15000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#pop1 has less samples (50%)
sample3 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=10000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#pop1 has less samples (25%)
sample4 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=5000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

save(sample1, sample2, sample3, sample4, file= "/work/yc16575/mr.trans/result/sim_difference_samplesize_290321.rdata")
#--------------------------sample size--------------------------#



#--------------------------MAF--------------------------#
#MAF POP1 = 0.25; MAF POP2 = 0.5
freq1 <- sample1 
freq2 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=20000, nid2=20000, maf1=0.25, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#MAF POP1 = 0.1; MAF POP2 = 0.5
freq3 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=20000, nid2=20000, maf1=0.1, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#MAF POP1 = 0.05; MAF POP2 = 0.5
freq4 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=20000, nid2=20000, maf1=0.05, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

save(freq1, freq2, freq3, freq4, file= "/work/yc16575/mr.trans/result/sim_difference_MAF_290321.rdata")
#--------------------------MAF--------------------------#



#--------------------------GENETIC CORRELATION--------------------------#
#genetic correlation = 0.25
corr1 <- sample1 
corr2 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#genetic correlation = 0.50
corr3 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0.5, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#genetic correlation = 0.75
corr4 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0.75, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#genetic correlation = 1.0
corr5 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=1, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)
save(corr1, corr2, corr3, corr4, corr5, file= "/work/yc16575/mr.trans/result/sim_difference_corr_290321.rdata")
#--------------------------GENETIC CORRELATION--------------------------#


#--------------------------HERITABITILIES--------------------------#
#heritability pop1 = 0.1; pop2 = 0.9
herit1 <- sample1 
herit2 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.1, hsq2=0.9, rg=0, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#heritability pop1 = 0.1; pop2 = 0.4
herit3 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.1, hsq2=0.4, rg=0, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#heritability pop1 = 0.1; pop2 = 0.1
herit4 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.1, hsq2=0.1, rg=0, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#heritability pop1 = 0.9; pop2 = 0.9
herit5 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.9, hsq2=0.9, rg=0, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)
save(herit1, herit2, herit3, herit4, herit5, file= "/work/yc16575/mr.trans/result/sim_difference_herit_290321.rdata")
#--------------------------HERITABITILIES--------------------------#


#--------------------------EFFECT ESTIMATES--------------------------#
effect1 <- sample1 

#Effect estimate pop1 = 0.8; pop2 = 0.8
effect2 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.8, biv2=0.8, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#Effect estimate pop1 = 0.2; pop2 = 0.8
effect3 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.8, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#Effect estimate pop1 = 0.2; pop2 = 0.2
effect4 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=-0.2, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#Effect estimate pop1 = 0.2; pop2 = 0.8
effect5 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=-0.8, inst_pval1=1, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)
save(effect1, effect2, effect3, effect4, effect5, file= "/work/yc16575/mr.trans/result/sim_difference_effect_290321.rdata")

#--------------------------EFFECT ESTIMATES--------------------------#


#--------------------------P VALUES--------------------------#
pval1 <- sample1 

#P value pop1 = 0.05; pop2 = 1
pval2 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=0.05, inst_pval2=1)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#P value pop1 = 0.05; pop2 = 0.05
pval3 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=0.05, inst_pval2=0.05)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#P value pop1 = 0.05/100; pop2 = 0.05
pval4 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=0.05/100, inst_pval2=0.05)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)

#P value pop1 = 0.05/100; pop2 = 0.05/100
pval5 <- mclapply(1:1000, function(i)
          {
            set.seed(i)
            res <- sim1(nsnp=100, hsq1=0.4, hsq2=0.4, rg=0, nid1=20000, nid2=20000, maf1=0.5, maf2=0.5, biv1=0.2, biv2=0.2, inst_pval1=0.05/100, inst_pval2=0.05/100)
            invisible(capture.output(s1 <- summary(res$model1, fit.measures=TRUE)))
            invisible(capture.output(s2 <- summary(res$model2, fit.measures=TRUE)))
            mod1 <- s1$FIT['aic']
            mod2 <- s2$FIT['aic']
            aic <- c(mod1=s1$FIT['aic'], mod2=s2$FIT['aic'])
            aic_diff <- mod1 - mod2
            return(list(s1=s1, s2=s2, aic=aic, aic_diff=aic_diff))
          }, mc.cores=16)
save(pval1, pval2, pval3, pval4, pval5, file= "/work/yc16575/mr.trans/result/sim_difference_pval_290321.rdata")


