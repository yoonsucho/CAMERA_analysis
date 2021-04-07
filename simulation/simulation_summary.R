
################## SUMMARY FOR THE SIMULATION ###################
#--------------------------sample size--------------------------#

load("/work/yc16575/mr.trans/result/sim_difference_samplesize_290321.rdata")

simres <- lapply(1:1000, function(x){
	sim <- x
	a <- sample1[[x]]$aic_diff
	b <- sample2[[x]]$aic_diff
	c <- sample3[[x]]$aic_diff
	d <- sample4[[x]]$aic_diff
	res <- tibble(sim, a, b, c, d)
	return(res)
}) %>% bind_rows 

png("/work/yc16575/mr.trans/result/sample1.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$a)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/sample2.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$b)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/sample3.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$c)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/sample4.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$d)
	abline(v=0)
dev.off()


#--------------------------MAF--------------------------#
rm(list=ls())
load("/work/yc16575/mr.trans/result/sim_difference_MAF_290321.rdata")

simres <- lapply(1:1000, function(x){
	sim <- x
	a <- freq1[[x]]$aic_diff
	b <- freq2[[x]]$aic_diff
	c <- freq3[[x]]$aic_diff
	d <- freq4[[x]]$aic_diff
	res <- tibble(sim, a, b, c, d)
	return(res)
}) %>% bind_rows 

png("/work/yc16575/mr.trans/result/freq1.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$a)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/freq2.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$b)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/freq3.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$c)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/freq4.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$d)
	abline(v=0)
dev.off()


#--------------------------GENETIC CORRELATION--------------------------#
rm(list=ls())
load("/work/yc16575/mr.trans/result/sim_difference_corr_290321.rdata")

simres <- lapply(1:1000, function(x){
	sim <- x
	a <- corr1[[x]]$aic_diff
	b <- corr2[[x]]$aic_diff
	c <- corr3[[x]]$aic_diff
	d <- corr4[[x]]$aic_diff
	e <- corr5[[x]]$aic_diff
	res <- tibble(sim, a, b, c, d, e)
	return(res)
}) %>% bind_rows 

png("/work/yc16575/mr.trans/result/corr1.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$a)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/corr2.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$b)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/corr3.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$c)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/corr4.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$d)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/corr5.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$e)
	abline(v=0)
dev.off()


#--------------------------HERITABITILIES--------------------------#
rm(list=ls())
load("/work/yc16575/mr.trans/result/sim_difference_herit_290321.rdata")

simres <- lapply(1:1000, function(x){
	sim <- x
	a <- herit1[[x]]$aic_diff
	b <- herit2[[x]]$aic_diff
	c <- herit3[[x]]$aic_diff
	res <- tibble(sim, a, b, c)
	return(res)
}) %>% bind_rows 

png("/work/yc16575/mr.trans/result/herit1.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$a)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/herit2.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$b)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/herit3.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$c)
	abline(v=0)
dev.off()




#--------------------------P VALUES--------------------------#
rm(list=ls())
load("/work/yc16575/mr.trans/result/sim_difference_pval_290321.rdata")

simres <- lapply(1:1000, function(x){
	sim <- x
	a <- pval1[[x]]$aic_diff
	b <- pval2[[x]]$aic_diff
	c <- pval3[[x]]$aic_diff
	d <- pval4[[x]]$aic_diff
	e <- pval5[[x]]$aic_diff
	res <- tibble(sim, a, b, c, d, e)
	return(res)
}) %>% bind_rows 

png("/work/yc16575/mr.trans/result/pval1.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$a)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/pval2.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$b)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/pval3.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$c)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/pval4.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$d)
	abline(v=0)
dev.off()

png("/work/yc16575/mr.trans/result/pval5.png", 
		width = 15, height = 12, units = "in", res = 300)
	hist(simres$e)
	abline(v=0)
dev.off()

