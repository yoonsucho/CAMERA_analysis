library(CAMeRa)
library(dplyr)
library(ggplot2)
library(ieugwasr)
select_api("private")


bfile_dir <- "/Users/gh13047/repo/opengwas-api-internal/opengwas-api/app/ld_files"
x <- CAMERA$new(
  exposure_ids=c(
    "ukb-e-23104_CSA", 
    "ukb-e-21001_AFR", 
    "ukb-b-19953", 
    "bbj-a-1"
  ), 
  outcome_ids=c(
    "ukb-e-411_CSA", 
    "ukb-e-411_AFR", 
    "ieu-a-7", 
    "bbj-a-109"
  ), 
  pops = c("SAS", "AFR", "EUR", "EAS"),
  bfiles=file.path(bfile_dir, c("SAS", "AFR", "EUR", "EAS")),
  plink = "plink",
  radius=50000, 
  clump_pop="EUR"
)

x$extract_instruments()
x$make_outcome_data()
x$extract_instrument_regions()
x$harmonise()
x$cross_estimate()

x <- readRDS(system.file("extdata", "example-CAMERA.rds", package="CAMeRa"))
x$exposure_ids
x$outcome_ids

x$instrument_fema
x$harmonise()
x$cross_estimate() %>% mutate(outcome="CHD") %>% {
    ggplot(., aes(x=Estimate, y=pops)) +
    geom_point(data=. %>% filter(pops != "All"), aes(colour=pops, size=Qjpval < 0.05)) +
    geom_errorbarh(data=. %>% filter(pops != "All"), aes(colour=pops, xmin=Estimate - 1.96 * `Std. Error`, xmax=Estimate + 1.96 * `Std. Error`), height=0) +
    geom_point(data=. %>% filter(pops == "All")) +
    geom_errorbarh(data=. %>% filter(pops == "All"), aes(xmin=Estimate - 1.96 * `Std. Error`, xmax=Estimate + 1.96 * `Std. Error`), height=0) +
    facet_grid(outcome ~ .) +
    scale_colour_brewer(type="qual", palette=3) +
    geom_vline(xintercept=0, linetype="dotted") +
    scale_y_discrete(limits=c(unique(.$pops)[unique(.$pops) != "All"], "All")) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    labs(y="", x="MR estimate of LDL cholesterol on stroke", colour="Ancestry", size="Heterogeneity\n-log10(p-value)")
}


x$estimate_instrument_heterogeneity_per_variant()
x$mrgxe()
x$mrgxe_plot()


CAMeRa::fixed_effects_meta_analysis(x$mrgxe_res$a, x$mrgxe_res$a_se) %>%
    {tibble(a=.$beta, se=.$se, pval=.$pval, nsnp=.$Qdf+1, Qpval=.$Qpval)}


x$pleiotropy()
temp <- x$pleiotropy_agreement %>% mutate(outcome="CHD") %>% bind_rows() %>%
    filter(metric=="Sign")
inner_join(subset(temp, datum=="Expected"), subset(temp, datum == "Observed"), by=c("disc", "rep", "outcome")) %>%
    ggplot(., aes(x=value.x, y=value.y)) +
    geom_point(aes(colour=paste(rep), shape=outcome)) +
    geom_smooth(method="lm", aes(colour=paste(disc)), se=FALSE) +
    facet_wrap(~ paste(disc)) +
    geom_abline(slope=1, intercept=0, linetype="dotted") +
    labs(x="Expected outlier replication", y="Observed outlier replication", colour="Ancestry", shape="Outcome") +
    scale_colour_brewer(type="qual", palette=3)


