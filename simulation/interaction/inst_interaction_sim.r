## Objective - The interaction model aims to maximise power when effects are shared, and estimate the heterogeneity when effects are distinct. Perform nested interaction model.

library(dplyr)
library(simulateGP)
library(MASS)
library(metafor)
library(lmtest)
library(TwoSampleMR)
library(parallel)

sim1 <- function(nsnp, hsq1, hsq2, rg, nid1, nid2, maf1, maf2, biv1, biv2, inst_pval1, inst_pval2)
{
  # Generate effects for exposures in pop1 and pop2
  s1 <- sqrt(hsq1 / nsnp)
  s2 <- sqrt(hsq2 / nsnp)
  sigma <- matrix(c(1, rg, rg, 1), 2, 2)
  eff <- MASS::mvrnorm(nsnp, c(0,0), sigma)
  
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
  
  index <- union(
    subset(ss1, pval.exposure < inst_pval1)$SNP,
    subset(ss2, pval.exposure < inst_pval2)$SNP
  )
  dat <- dat[index,]
  return(dat)
}

make_geno <- function(nid, nsnp, af)
{
    if(length(af) == 1)
    {
        return(matrix(rbinom(nid * nsnp, 2, af), nid, nsnp))
    } else if(length(af) == nsnp)
    {
        do.call(cbind, lapply(af, function(p) rbinom(nid, 2, p))) %>% return
    } else {
        stop("af should be length 1 or nsnp")
    }
}

sim1 <- function(nsnp, hsq, sigma, nid, af, biv, inst_pval) {
    # input checks
    npop <- length(hsq)
    stopifnot(npop == length(nid))
    stopifnot(npop == length(af))
    stopifnot(npop == length(biv))
    stopifnot(npop == length(inst_pval))
    stopifnot(nrow(sigma) == ncol(sigma))
    stopifnot(nrow(sigma) == npop)

    map <- lapply(1:npop, function(i) tibble(snp=1:nsnp, af = af[[i]]))

    # Generate effects for exposures in pop1 and pop2
    eff <- MASS::mvrnorm(nsnp, rep(0, npop), sigma)

    # create genotypes
    g <- lapply(1:npop, function(i) make_geno(nid[i], nsnp, af[[i]]))

    # create confounder
    u <- lapply(1:npop, function(i) rnorm(nid[i], sd=sqrt(0.1)))
  
    # prs
    prs <- lapply(1:npop, function(i) scale(g[[i]] %*% eff[,i]) * sqrt(hsq[i]))

    # e
    e <- lapply(1:npop, function(i) rnorm(nid[i], 0, sqrt(1 - var(prs[[i]]) - var(u[[i]]))) )

    # create x
    x <- lapply(1:npop, function(i) prs[[i]] + u[[i]] + e[[i]] )

    # create y
    y <- lapply(1:npop, function(i) x[[i]] * biv[i] + u[[i]] + rnorm(nid[i], 0, sqrt(1 - biv[i]^2 * var(x[[i]]) - var(u[[i]]))))

    # Summary data
    ss <- lapply(1:npop, function(i) get_effs(x[[i]], y[[i]], g[[i]]))

    # Get SNPs to include
    index <- Reduce(union, lapply(1:npop, function(i) subset(ss[[i]], pval.exposure < inst_pval[i])$SNP))
  
    ss <- lapply(ss, function(x) x[index, ])
    return(ss)
}

int_analysis <- function(dat){
    npop <- length(dat)
    d <- lapply(1:npop, function(i) tibble(snp=dat[[i]]$SNP, x=dat[[i]]$beta.exposure, y=dat[[i]]$beta.outcome, xse=dat[[i]]$se.exposure, yse=dat[[i]]$se.outcome, pop=i)) %>% bind_rows()
    mod1 <- lm(y ~ -1 + x, data=d, weight=1/d$yse^2)
    mod2 <- lm(y ~ -1 + x:as.factor(pop), data=d, weight=1/d$yse^2)
    smod1 <- summary(mod1)
    smod2 <- summary(mod2)
    modcomp <- anova(mod1, mod2)$P[2]

    mod3 <- lm(y ~ -1 + x + x:as.factor(pop), data=d, weight=1/d$yse^2)
    smod3 <- summary(mod3)
    modcomp2 <- anova(mod1, mod3)$P[2]

    jt <- lmtest::jtest(mod1, mod2)$P[1]
    ct <- lmtest::coxtest(mod1, mod2)$P[1]
    aicdif <- AIC(mod2)-AIC(mod1)

    o <- lapply(dat, function(d) mr(d, method=c("mr_wald_ratio", "mr_ivw"))) %>% bind_rows()
    r <- metafor::rma.uni(o$b, sei=o$se, method="FE")

    pop_ests <- bind_rows(
        tibble(
            method="interaction",
            pop=1:npop,
            b=smod2$coef[,1],
            se=smod2$coef[,2],
            pval=smod2$coef[,4]
        ),
        tibble(
            method="ivw",
            pop=1:npop,
            b=o$b,
            se=o$se,
            pval=o$pval,
            nsnp=o$nsnp
        )
    ) %>% mutate(pop=as.character(pop))

    bind_rows(
        tibble(
            method="nested anova",
            pop="all",
            b=smod1$coef[1,1],
            se=smod1$coef[1,2],
            pval=smod1$coef[1,4],
            diff_pval=modcomp2,
            npop=npop
        ),
        tibble(
            method="non-nested anova",
            pop="all",
            b=smod1$coef[1,1],
            se=smod1$coef[1,2],
            pval=smod1$coef[1,4],
            diff_pval=modcomp,
            npop=npop
        ),
        tibble(
            method="jtest",
            pop="all",
            b=smod1$coef[1,1],
            se=smod1$coef[1,2],
            pval=smod1$coef[1,4],
            diff_pval=jt,
            npop=npop
        ),
        tibble(
            method="coxtest",
            pop="all",
            b=smod1$coef[1,1],
            se=smod1$coef[1,2],
            pval=smod1$coef[1,4],
            diff_pval=ct,
            npop=npop
        ),
        tibble(
            method="AIC",
            pop="all",
            b=smod1$coef[1,1],
            se=smod1$coef[1,2],
            pval=smod1$coef[1,4],
            diff_pval=aicdif,
            npop=npop
        ),
        tibble(
            method="FE meta analysis",
            pop="all",
            b=r$beta[1],
            se=r$se,
            pval=r$pval,
            diff_pval=r$QEp,
            npop=length(dat)
        )
    ) %>% bind_rows(pop_ests) %>% return

}

# a <- sim1(nsnp=100, hsq=c(0.3, 0.4), sigma=matrix(c(1,0.2, 0.2, 1), 2, 2), nid=c(20000, 10000), af=c(0.2, 0.3), c(0.1, 0.12), c(1,1))
# int_analysis(a)
# sigma <- matrix(1, 3, 3)
# diag(sigma) <- 1
# b <- sim1(nsnp=100, hsq=c(0.3, 0.4, 0.3), sigma=sigma, nid=c(20000, 10000, 20000), af=c(0.2, 0.3, 0.3), c(0.1, 0.12, -0.1), c(5e-8, 5e-8, 5e-8))
# int_analysis(b)


# Analytical procedure
# - Overall estimate
# - Per ancestry estimate
# - ANOVA to determine if per ancestry offers better fit

# Scenarios

# - Power is lower in one population
# - FDR under null

param <- expand.grid(
    nsnp=100,
    npop=c(2,3,4),
    hsq=0.1,
    rg=c(0, 0.5, 1),
    biv_m=c(0, 1),
    biv_sd=c(0, 0.1, 0.2),
    sim=1:100
)
param$index <- 1:nrow(param)

res <- mclapply(1:nrow(param), function(i)
{
    sigma <- matrix(param$rg[i], param$npop[i], param$npop[i])
    diag(sigma) <- 1
    af <- lapply(1:param$npop[i], function(i) runif(param$nsnp[i], 0.01, 0.99))
    sim1(nsnp=param$nsnp[i], hsq=rep(param$hsq[i], param$npop[i]), sigma=sigma, nid=rep(10000, param$npop[i]), af=af, biv=rnorm(param$npop[i], param$biv_m[i], param$biv_sd[i]), rep(5e-8, param$npop[i])) %>%
    int_analysis %>% dplyr::select(-c(npop, nsnp)) %>% bind_cols(param[i,], .)
}, mc.cores=40)
save(res, file="inst_interaction_sim1.rdata")

