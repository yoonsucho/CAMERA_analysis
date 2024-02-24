library(here)
library(CAMeRa)
library(dplyr)
library(ieugwasr)
select_api("private")

mntdir <- "/Volumes/010/working/data/examples"
fn <- tibble(
    fp = list.files(file.path(mntdir, "results"), full.names=TRUE),
    eur_id_x = NA,
    eas_id_x = NA,
    eur_id_y = NA,
    eas_id_y = NA
)


for(i in 1:nrow(fn)) {
    message(i)
    load(fn$fp[i])
    fn$eur_id_x[i] <- a$exposure_ids[1]
    fn$eas_id_x[i] <- a$exposure_ids[2]
    fn$eur_id_y[i] <- a$outcome_ids[1]
    fn$eas_id_y[i] <- a$outcome_ids[2]
    x <- CAMERA$new()
    x$import(a)
    for(j in 1:length(x$instrument_regions)) {
        if(!is.null(x$instrument_regions[[j]]))
        {
            x$instrument_regions[[j]][[1]] <- CAMeRa:::generate_vid(x$instrument_regions[[j]][[1]])
            x$instrument_regions[[j]][[2]] <- CAMeRa:::generate_vid(x$instrument_regions[[j]][[2]])
        }
    }
    x$fema_regional_instruments()
    x$make_outcome_data(exp=x$instrument_fema)
    x$instrument_raw <- CAMeRa:::generate_vid(x$instrument_raw)
    x$make_outcome_data(exp=x$instrument_raw)
    x$source <- "OpenGWAS"
    x$harmonise(exp=x$instrument_fema)
    x$cross_estimate()
    # x$estimate_instrument_heterogeneity_per_variant()
    # x$mrgxe()
    # x$pleiotropy()
    save(x, file=here("data", "bbj", basename(fn$fp[i])))
}

fn$newpath <- file.path(here("data", "bbj"), basename(fn$fp))
saveRDS(fn, here("data", "bbj", "fn.rds"))

i <- 1
res <- list()
for(i in 1:nrow(fn)) {
    r <- list()
    message(i)
    load(fn$newpath[i])
    x$harmonise(exp=x$instrument_raw)
    o1 <- x$cross_estimate() %>% mutate(f=basename(fn$fp[i]), instrument="Raw")
    x$harmonise(exp=x$instrument_fema)
    o2 <- x$cross_estimate() %>% mutate(f=basename(fn$fp[i]), instrument="FEMA")
    r$mrres <- bind_rows(o1, o2)

    o1 <- x$estimate_instrument_specificity(instrument=x$instrument_raw) %>% mutate(instrument="Raw")
    o2 <- x$estimate_instrument_specificity(instrument=x$instrument_fema) %>% mutate(instrument="FEMA")
    r$instrument_specificity <- bind_rows(o1, o2)
    p1 <- x$instrument_heterogeneity(instrument=x$instrument_raw) %>% mutate(instrument="Raw")
    p2 <- x$instrument_heterogeneity(instrument=x$instrument_fema) %>% mutate(instrument="FEMA")
    r$instrument_heterogeneity <- bind_rows(p1, p2)

    x$estimate_instrument_heterogeneity_per_variant()
    x$mrgxe()
    r$mrgxe <- CAMeRa::fixed_effects_meta_analysis(x$mrgxe_res$a, x$mrgxe_res$a_se) %>%
    {tibble(a=.$beta, se=.$se, pval=.$pval, nsnp=.$Qdf+1, Qpval=.$Qpval)}

    x$pleiotropy()
    x$pleiotropy_outliers %>% str
    x$pleiotropy_Q_outliers %>% filter(p.adjust(Qpval, "fdr") < 0.05)
    r$pleiotropy <- x$pleiotropy_agreement
    res[[i]] <- r
}

saveRDS(res, file=here("data", "bbj", "res.rds"))
