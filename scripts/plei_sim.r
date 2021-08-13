library(simulateGP)
library(tidyverse)
library(lavaan)
library(RadialMR)

sim1 <- function(nsnp, nid1x, nid2x, nid1y, nid2y, h2, S, biv1, biv2, pl1, pl2, absx)
{
	# Create map of variants
	map1 <- tibble(snp=1:nsnp, af=runif(nsnp))
	map2 <- tibble(snp=1:nsnp, af=runif(nsnp))

	# Create SNP-exposure effects
	x1 <- generate_gwas_params(map1, h2=h2, Pi=1, S=S)
	if(absx)
	{
		x1$beta <- abs(x1$beta)
	}
	x2 <- map2
	x2$beta <- x1$beta

	# Create SNP-outcome effects
	y1 <- x1
	y1$beta <- y1$beta * biv1 + pl1
	y2 <- x2
	y2$beta <- y2$beta * biv2 + pl2

	x1s <- generate_gwas_ss(x1, nid=nid1x)
	x2s <- generate_gwas_ss(x2, nid=nid2x)
	y1s <- generate_gwas_ss(y1, nid=nid1y)
	y2s <- generate_gwas_ss(y2, nid=nid2y)

	dat <- tibble(
		y1 = y1s$bhat,
		y2 = y2s$bhat,
		x1 = x1s$bhat,
		x2 = x2s$bhat,
		yse1 = y1s$se,
		yse2 = y2s$se,
		xse1 = x1s$se,
		xse2 = x2s$se,
		r1 = y1/x1,
		r2 = y2/x2,
		w1 = sqrt(x1^2 / yse1^2),
		w2 = sqrt(x2^2 / yse2^2),
		o1 = r1 * w1,
		o2 = r2 * w2
	)

	out <- list()
	out$ivw1 <- TwoSampleMR::mr_ivw(dat$x1, dat$y1, dat$xse1, dat$yse1)$b
	out$ivw2 <- TwoSampleMR::mr_ivw(dat$x2, dat$y2, dat$xse2, dat$yse2)$b
	out$rm1 <- summary(lm(o1 ~ -1 + w1, data=dat))$coef[1,1]
	out$rm2 <- summary(lm(o2 ~ -1 + w2, data=dat))$coef[1,1]

	model1 <- '
	y1 ~ biv*x1
	y2 ~ biv*x2
	'
	mod1 <- sem(model1, data=dat)
	invisible(capture.output(mod1 <- summary(mod1)))
	out$sem11 <- mod1$PE$est[1]
	out$sem12 <- mod1$PE$est[2]

	model2 <- '
	y1 ~ biv_1*x1
	y2 ~ biv_2*x2
	'
	mod2 <- sem(model2, data=dat)
	invisible(capture.output(mod2 <- summary(mod2)))
	out$sem21 <- mod2$PE$est[1]
	out$sem22 <- mod2$PE$est[2]

	model3 <- '
	o1 ~ biv*w1
	o2 ~ biv*w2
	'
	mod3 <- sem(model3, data=dat)
	invisible(capture.output(mod3 <- summary(mod3)))
	out$sem31 <- mod3$PE$est[1]
	out$sem32 <- mod3$PE$est[2]
	model4 <- '
	o1 ~ biv_1*w1
	o2 ~ biv_2*w2
	'
	mod4 <- sem(model4, data=dat)
	invisible(capture.output(mod4 <- summary(mod4)))
	out$sem41 <- mod4$PE$est[1]
	out$sem42 <- mod4$PE$est[2]
	return(as_tibble(out))
}



