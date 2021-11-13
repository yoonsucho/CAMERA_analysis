library(simulateGP)
library(TwoSampleMR)

nsnp <- 50
nid <- 100000
g <- make_geno(nid, nsnp, af=0.3)
b1 <- choose_effects(nsnp, 0.3)
b2 <- choose_effects(nsnp, 0.3)

bxy1 <- 0.1
bxy2 <- 0.1

x1 <- make_phen(b1, g)
x2 <- make_phen(b1, g)
y1 <- make_phen(c(bxy1), cbind(x1))
y2 <- make_phen(c(bxy2), cbind(x2))




x <- get_effs(x1, x2, g)
y <- get_effs(y1, y2, g)

mr(x, method="mr_ivw")
mr(y, method="mr_ivw")




x2a <- x2 * 6
y2a <- y2 * 4

y1a <- y1 * 10

mr(get_effs(y1a, y2a, g), meth="mr_ivw")$b 

mr(get_effs(x1, y1, g), meth="mr_ivw")
mr(subset(get_effs(x2, y2, g), pval.exposure < 5e-8), meth="mr_ivw")

mr(subset(get_effs(x2, y2a, g), pval.exposure < 5e-8), meth="mr_ivw")

mr(subset(get_effs(x1, y1, g), pval.exposure < 5e-8), meth="mr_ivw")

mr(subset(get_effs(y1, y2a, g)), meth="mr_ivw")


library(tidyverse)
nsim <- 100
p <- tibble(
	xvar1 = runif(nsim, 0.1, 10),
	xvar2 = runif(nsim, 0.1, 10),
	yvar1 = runif(nsim, 0.1, 10),
	yvar2 = runif(nsim, 0.1, 10),
	bxy1 = runif(nsim, 0.2, 0.6) * sample(c(1, -1), nsim, rep=T),
	bxy2 = runif(nsim, 0.2, 0.6) * sample(c(1, -1), nsim, rep=T),
	bxy1hat = NA,
	bxy2hat = NA,
	xscale = NA,
	yscale = NA
)

nsnp <- 50
nid <- 100000
g1 <- make_geno(nid, nsnp, af=0.3)
g2 <- make_geno(nid, nsnp, af=0.3)
b1 <- choose_effects(nsnp, 0.3)
b2 <- choose_effects(nsnp, 0.3)

for(i in 1:nsim)
{
	message(i)
	x1 <- make_phen(b1, g1)
	x2 <- make_phen(b1, g1)
	y1 <- make_phen(c(p$bxy1[i], b2), cbind(x1, g2))
	y2 <- make_phen(c(p$bxy2[i], b2), cbind(x2, g2))

	x1a <- x1 * p$xvar1[i]
	x2a <- x2 * p$xvar2[i]
	y1a <- y1 * p$yvar1[i]
	y2a <- y2 * p$yvar2[i]

	p$xscale[i] <- get_effs(x1a, x2a, g1) %>%
		subset(., pval.exposure < 5e-8) %>%
		mr(., meth="mr_ivw") %>%
		{.$b}
	p$yscale[i] <- get_effs(y1a, y2a, g2) %>%
		subset(., pval.exposure < 5e-8) %>%
		mr(., meth="mr_ivw") %>%
		{.$b}

	p$bxy1hat[i] <- get_effs(x1a, y1a, g1) %>%
		subset(., pval.exposure < 5e-8) %>%
		mr(., meth="mr_ivw") %>%
		{.$b}

	p$bxy2hat[i] <- get_effs(x2a, y2a, g1) %>%
		subset(., pval.exposure < 5e-8) %>%
		mr(., meth="mr_ivw") %>%
		{.$b}
}




nsim <- 100
p2 <- tibble(
	xvar1 = runif(nsim, 0.1, 10),
	xvar2 = runif(nsim, 0.1, 10),
	yvar1 = runif(nsim, 0.1, 10),
	yvar2 = runif(nsim, 0.1, 10),
	bxy1 = runif(nsim, 0.2, 0.6) * sample(c(1, -1), nsim, rep=T),
	bxy2 = bxy1,
	bxy1hat = NA,
	bxy2hat = NA,
	xscale = NA,
	yscale = NA
)

for(i in 1:nsim)
{
	message(i)
	x1 <- make_phen(b1, g1)
	x2 <- make_phen(b1, g1)
	y1 <- make_phen(c(p2$bxy1[i], b2), cbind(x1, g2))
	y2 <- make_phen(c(p2$bxy2[i], b2), cbind(x2, g2))

	x1a <- x1 * p2$xvar1[i]
	x2a <- x2 * p2$xvar2[i]
	y1a <- y1 * p2$yvar1[i]
	y2a <- y2 * p2$yvar2[i]

	p2$xscale[i] <- get_effs(x1a, x2a, g1) %>%
		subset(., pval.exposure < 5e-8) %>%
		mr(., meth="mr_ivw") %>%
		{.$b}
	p2$yscale[i] <- get_effs(y1a, y2a, g2) %>%
		subset(., pval.exposure < 5e-8) %>%
		mr(., meth="mr_ivw") %>%
		{.$b}

	p2$bxy1hat[i] <- get_effs(x1a, y1a, g1) %>%
		subset(., pval.exposure < 5e-8) %>%
		mr(., meth="mr_ivw") %>%
		{.$b}

	p2$bxy2hat[i] <- get_effs(x2a, y2a, g1) %>%
		subset(., pval.exposure < 5e-8) %>%
		mr(., meth="mr_ivw") %>%
		{.$b}
}







ggplot(p, aes(x=bxy1/bxy2, y=bxy1hat/bxy2hat)) +
geom_point()

ggplot(p, aes(x=bxy1/bxy2, y=bxy1hat/bxy2hat * yscale/xscale)) +
geom_point()

ggplot(p2, aes(x=bxy1, y=bxy1hat)) +
geom_point()


ggplot(p, aes(x=bxy2, y=bxy2hat * xvar2/yvar2)) +
geom_point()

