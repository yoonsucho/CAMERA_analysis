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


# ieu-b-25 ieu-b-5… ieu-b-38 bbj-a-52
x <- CAMERA$new(
  exposure_ids=c(
    "ieu-b-25", 
    "ieu-b-5071"
  ), 
  outcome_ids=c(
    "ukb-b-20175", 
    "bbj-a-52"
  ), 
  pops = c("EUR", "EAS"),
  radius=50000, 
  clump_pop="EUR"
)
x$extract_instruments()
x$make_outcome_data()
x$harmonise()
x$cross_estimate()
x$extract_instrument_regions()
x$fema_regional_instruments()
x$make_outcome_data()
x$harmonise(exp=x$instrument_fema)
x$cross_estimate()

fn2 <- tibble(
    eur_id_x = "ieu-b-25",
    eas_id_x = "ieu-b-5071",
    eur_id_y = "ukb-b-20175",
    eas_id_y = "bbj-a-52",
    newpath = here("data", "bbj", "examples_cig_sbp_new.rdata")
)
save(x, file=here("data", "bbj", basename(fn2$newpath)))
fn <- bind_rows(fn, fn2)


x <- CAMERA$new(
  exposure_ids=c(
    "ieu-b-25", 
    "ieu-b-5071"
  ), 
  outcome_ids=c(
    "ukb-b-20175", 
    "ieu-b-5075"
  ), 
  pops = c("EUR", "EAS"),
  radius=50000, 
  clump_pop="EUR"
)
x$extract_instruments()
x$make_outcome_data()
x$harmonise()
x$cross_estimate()
x$extract_instrument_regions()
x$fema_regional_instruments()
x$make_outcome_data()
x$harmonise(exp=x$instrument_fema)
x$cross_estimate()

fn2 <- tibble(
    eur_id_x = "ieu-b-25",
    eas_id_x = "ieu-b-5071",
    eur_id_y = "ukb-b-20175",
    eas_id_y = "ieu-b-5075",
    newpath = here("data", "bbj", "examples_cig_sbp_new2.rdata")
)
save(x, file=here("data", "bbj", basename(fn2$newpath)))
fn <- bind_rows(fn, fn2)



i <- 1
res <- list()
for(i in 60:nrow(fn)) {
    r <- list()
    message(i)
    load(fn$newpath[i])
    x$harmonise(exp=x$instrument_raw)
    o1 <- x$cross_estimate() %>% mutate(f=basename(fn$newpath[i]), instrument="Raw")
    x$harmonise(exp=x$instrument_fema)
    o2 <- x$cross_estimate() %>% mutate(f=basename(fn$newpath[i]), instrument="FEMA")
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


gi <- gwasinfo(unique(c(fn$eur_id_x, fn$eur_id_y)))
gi$trait

.simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}
gi$trait <- sapply(gi$trait, .simpleCap)
gi$trait[grepl("Breast Cancer", gi$trait)] <- "Breast Cancer"
gi$trait[grepl("Systolic", gi$trait)] <- "Systolic Blood Pressure"
gi$trait

fn <- left_join(fn, dplyr::select(gi, id, exposure=trait), by=c("eur_id_x"="id")) %>%
    left_join(., dplyr::select(gi, id, outcome=trait), by=c("eur_id_y"="id"))

saveRDS(res, file=here("data", "bbj", "res.rds"))
saveRDS(fn, here("data", "bbj", "fn.rds"))

