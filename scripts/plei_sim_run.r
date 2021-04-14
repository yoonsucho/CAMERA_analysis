library(simulateGP)
library(tidyverse)
library(lavaan)
library(RadialMR)
library(parallel)
library(stringr)
source("plei_sim.r")

param <- expand.grid(
	nsnp = 100,
	nid1x = 100000,
	nid2x = 100000,
	nid1y = 100000,
	nid2y = 100000,
	h2 = 0.5,
	S = 1,
	biv1 = seq(-0.2, 0.2, by=0.02),
	biv2 = seq(-0.2, 0.2, by=0.02),
	pl1 = c(0, 0.2),
	pl2 = c(0, 0.2),
	absx = c(FALSE),
	sim=1
) %>% as_tibble()
dim(param)
res <- mclapply(1:nrow(param), function(i)
{
	message(i)
	do.call(sim1, args=param[i,] %>% select(-c(sim))) %>%
		bind_cols(param[i, ], .)
}, mc.cores=16) %>% bind_rows()

p1 <- res %>%
	pivot_longer(cols=c(ivw1,ivw2,rm1,rm2,sem11,sem12,sem21,sem22,sem31,sem32,sem41,sem42), names_to="method", values_to="bhat") %>%
	pivot_longer(cols=c(biv1, biv2), names_to="pop", values_to="biv") %>%
	pivot_longer(cols=c(pl1, pl2), names_to="pop2", values_to="pl") %>%
	mutate(pop=gsub("biv", "", pop), pop2=gsub("pl", "", pop), methodest=str_sub(method, -1)) %>% 
	dplyr::filter(pop==pop2 & pop == methodest) %>%
	dplyr::select(-c(pop2, methodest)) %>%
	group_by(pop, biv, pl, absx, method) %>%
	summarise(
		se=sd(bhat, na.rm=T),
		bhat=mean(bhat, na.rm=T),
		n=n()
	) %>%
	ggplot(., aes(x=biv, y=bhat)) +
	geom_point(aes(colour=pop, size=abs(bhat-biv))) +
	geom_errorbar(width=0, aes(ymin=bhat-1.96*se, ymax=bhat+1.96*se, colour=pop)) +
	facet_grid(method ~ pl)
ggsave(p1, file="temp.pdf")


p2 <- res %>%
	pivot_longer(cols=c(ivw1,ivw2,rm1,rm2,sem11,sem12,sem21,sem22,sem31,sem32,sem41,sem42), names_to="method", values_to="bhat") %>%
	group_by(biv1, biv2, pl1, pl2, absx, method) %>%
	summarise(
		se=sd(bhat, na.rm=T),
		bhat=mean(bhat, na.rm=T),
		n=n()
	) %>%
	ggplot(., aes(x=biv1, y=biv2)) +
	geom_line(aes(colour=))
	geom_point(aes(alpha=abs(bhat-biv2))) +
	facet_grid(method ~ pl1*pl2)
ggsave(p2, file="temp2.pdf")






set.seed(123)
nsnp <- 100

# Create map of variants
map1 <- tibble(snp=1:nsnp, af=runif(100))
map2 <- tibble(snp=1:nsnp, af=runif(100))

# Create SNP-exposure effects
x1 <- generate_gwas_params(map1, h2=0.5, Pi=1, S=0)
x2 <- x1
x2$beta <- x1$beta

# Create SNP-outcome effects
y1 <- x1
y1$beta <- y1$beta * 0.1
y2 <- x2
y2$beta <- y2$beta * 0.1

x1s <- generate_gwas_ss(x1, nid=100000)
x2s <- generate_gwas_ss(x2, nid=100000)
y1s <- generate_gwas_ss(y1, nid=100000)
y2s <- generate_gwas_ss(y2, nid=100000)

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


TwoSampleMR::mr_ivw(dat$x1, dat$y1, dat$xse1, dat$yse1)
TwoSampleMR::mr_ivw(dat$x2, dat$y2, dat$xse2, dat$yse2)

ivw_radial(format_radial(dat$x1, dat$y1, dat$xse1, dat$yse1), weights=1)
ivw_radial(format_radial(dat$x2, dat$y2, dat$xse2, dat$yse2), weights=1)


model1 <- '
y1 ~ biv*x1
y2 ~ biv*x2
'
mod1 <- sem(model1, data=dat)
summary(mod1)

model2 <- '
y1 ~ biv_1*x1
y2 ~ biv_2*x2
'
mod2 <- sem(model2, data=dat)
summary(mod2)


# Radial MR model

model3 <- '
o1 ~ biv*w1
o2 ~ biv*w2
'
mod3 <- sem(model3, data=dat)
summary(mod3)


model4 <- '
o1 ~ biv_1*w1
o2 ~ biv_2*w2
'
mod4 <- sem(model4, data=dat)
summary(mod4)






### Introduce pleiotropy


set.seed(123)
nsnp <- 100

# Create map of variants
map1 <- tibble(snp=1:nsnp, af=runif(100))
map2 <- tibble(snp=1:nsnp, af=runif(100))

x_index <- 1:90
z_index <- 91:100

# Create SNP-exposure effects
x1 <- generate_gwas_params(map1, h2=0.5, Pi=1, S=0)
x1$beta <- abs(x1$beta)
x1$beta[z_index] <- 0
x2 <- x1
x2$beta <- x1$beta

# Create SNP-outcome effects

y1 <- x1
y1$beta <- y1$beta * 0.1 + 0.2
y2 <- x2
y2$beta <- y2$beta * 0.1 + 0.2

x1s <- generate_gwas_ss(x1, nid=100000)
x2s <- generate_gwas_ss(x2, nid=100000)
y1s <- generate_gwas_ss(y1, nid=100000)
y2s <- generate_gwas_ss(y2, nid=100000)

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

dat$x1[z_index] <- 0

plot(y1 ~ x1, dat, ylim=c(0, 0.3))
abline(lm(y1 ~ x1, dat))
abline(lm(y1 ~ -1 + x1, dat))

TwoSampleMR::mr_ivw(dat$x1, dat$y1, dat$xse1, dat$yse1)
TwoSampleMR::mr_ivw(dat$x2, dat$y2, dat$xse2, dat$yse2)
summary(lm(o1 ~ -1 + w1, data=dat))
summary(lm(o2 ~ -1 + w2, data=dat))


model1 <- '
y1 ~ biv*x1
y1 ~ plei*x2
'
mod1 <- sem(model1, data=dat)
summary(mod1)

model2 <- '
y1 ~ biv_1*x1
y2 ~ biv_2*x2
'
mod2 <- sem(model2, data=dat)
summary(mod2)





# Radial MR model

model3 <- '
o1 ~ biv*w1
o2 ~ biv*w2
'
mod3 <- sem(model3, data=dat)
summary(mod3)


model4 <- '
o1 ~ biv_1*w1
o2 ~ biv_2*w2
'
mod4 <- sem(model4, data=dat)
summary(mod4)










dat$r1 <- dat$y1/dat$x1
dat$r2 <- dat$y2/dat$x2
dat$w1 <- sqrt(dat$x1^2 / dat$yse1^2)
dat$w2 <- sqrt(dat$x2^2 / dat$yse2^2)
dat$o1 <- dat$r1 * dat$w1
dat$o2 <- dat$r2 * dat$w2


ivw_radial(format_radial(dat$x1, dat$y1, dat$xse1, dat$yse1), weights=1)



model3 <- '
o1 ~ biv_1*w1
o2 ~ biv_2*w2
'
mod3 <- sem(model3, data=dat)
summary(mod3)



