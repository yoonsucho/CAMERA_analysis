# Trans ancestral design for Mendelian Randomization
Develop trans-ancestral model for MR

## Background
Emerging availability of large scale GWAS in non-European populations has enabled `trans-ancestral' studies that are conducted using more than one population. Trans-ancestral study design have been suggested to overcome existing issues with limited sample sizes available for many diseases. 

There have been a few studies that introduced trans-ancestral design to MR. However, no guideline exist for conducting and reporting trans-ancestral MR so far. Here, we simulated possible scenarios for trans-ancestral MR to provide a guideline for the trans-ancestral MR study design.

### Can we perform two sample MR using the samples from different ancestry?
MR studies require ethnically homogenous samples. The following scenario can be developed:
- Size and direction of genetic effect can be different across the ancestries
- Causal structures are different across the ancestries.
- Minor allele frequencies are different
- Ancestry specific SNPs

### Simulation
1. MR of exposure in ASN against outcome in EUR or vice versa
-	Set different MAF, Effect size (e.g. How big the difference should be?)
- Suggested model: Y1,Y2 ~ X1 + X2, where X1 and Y1 is the exposure and the outcome in ASN, and X2 and Y2 is the expousre and the outcome in EUR, respectively.
2. MR of exposure and outcome in ASN using EUR discovery SNPs
3. Using SEM model to get an estimate of pleiotropy 
4. Using fine-mapping to help improve analyses

1. MR of exposure in ASN against outcome in EUR or vice versa
```{r}
#Sample (n = 2): Y1(N=1000), Y2(N=1000)
#Variants (n = 200): g1 … g200; MAF = 0.5
#100 have an effect on X in both pop 1 and pop 2
#50 only on pop 1 = X1
#50 only on pop 2 = X2
#X has the same causal effect on Y in both pop 1 and pop 2 = XB

#Need to generate genotype data for all 200 SNPs for two populations
#Need to generate X and Y for two populations
#Generate GWAS summary data of all 200 SNPs on X and Y for each population

#y1,y2 ~ x1 + x2 + xb

library(simulateGP)
effects1 <- effects2 <- choose_effects(200, 0.4)
effects1[1:50] <- 0
effects2[51:100] <- 0

# create genotypes
g1 <- make_geno(10000, 200, 0.5)
g2 <- make_geno(10000, 200, 0.5)

# create confounder
u1 <- rnorm(10000)
u2 <- rnorm(10000)

# create x
x1 <- make_phen(c(0.1, effects1), cbind(u1, g1))
x2 <- make_phen(c(0.1, effects2), cbind(u2, g2))

# create y
y1 <- make_phen(c(0.1, 0.2), cbind(u1, x1))
y2 <- make_phen(c(0.1, 0.2), cbind(u2, x2))

# Summary data for pop1 
ss1 <- get_effs(x1, y1, g1)

# Summary data for pop2
ss2 <- get_effs(x2, y2, g2)

# Organise the data
library(dplyr)
dat <- tibble(
	y1 = ss1$beta.outcome,
	y2 = ss2$beta.outcome,
	x1 = ss1$beta.exposure,
	x2 = ss2$beta.exposure
)

dat$xb <- rowMeans(cbind(dat$x1, dat$x2))
dat$xb[1:100] <- 0
dat$x1[1:50] <- 0
dat$x1[101:200] <- 0
dat$x2[51:100] <- 0
dat$x2[101:200] <- 0

# fit model
summary(lm(cbind(y1, y2) ~ xb, data=dat))
summary(lm(cbind(y1, y2) ~ x1, data=dat)); Any effect on y2 should be pleiotropy
summary(lm(cbind(y1, y2) ~ x2, data=dat))
summary(lm(cbind(y1, y2) ~ xb + x1 + x2, data=dat))
```



