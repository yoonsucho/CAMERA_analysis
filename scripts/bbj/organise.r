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


for(i in 25:nrow(fn)) {
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


