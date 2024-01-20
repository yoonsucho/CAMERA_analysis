a <- out$region_extract[[1]][[1]]
a
join


subset(a, trait == "LDL", select=c(vid, pop, beta)) %>% tidyr::pivot_wider(id_cols=vid, names_from=pop, values_from=beta)

make_matrices <- function(a) {
    traits <- unique(a$trait)
    pops <- unique(a$pop)
    lapply(traits, \(tr) {
        lapply(c("beta", "se", ))
        subset(a, trait == tr, select=c("vid")) %>%

    })

    for(i in 1:length(traits)) {
        l1[[i]] <- list()
        for(j in 1:length(pops)) {

        }
    }
}


z_meta_analysis <- function(beta_mat, se_mat, n, eaf_mat) {
    z_mat <- beta_mat / se_mat
    w_mat <- t(t(sqrt(eaf_mat * (1 - eaf_mat) * 2)) * sqrt(n))
    zw <- rowSums(z_mat * w_mat, na.rm=TRUE) / sqrt(rowSums(w_mat^2, na.rm=TRUE))
    pval <- qnorm(abs(zw), lower.tail = FALSE)
    return(pval)
}


