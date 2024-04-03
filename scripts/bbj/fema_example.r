library(dplyr)
library(here)
library(CAMeRa)
library(ggplot2)

load(here("data", "bbj", "examples_sbp_t2d.rdata"))
ls()

i <- 6
x$plot_regional_instruments(names(x$instrument_regions)[i]) + facet_grid(pop ~ region, scale="free_y")


p <- function (region, instrument_regions = x$instrument_regions, meta_analysis_regions = x$instrument_fema_regions) {
    r <- instrument_regions[[region]]
    d <- meta_analysis_regions[[region]]
    r <- lapply(r, function(y) y %>% dplyr::mutate(z = abs(beta)/se))
    r$fema <- d
    temp <- lapply(names(r), function(y) r[[y]] %>% dplyr::mutate(pop = y)) %>% 
        dplyr::bind_rows() %>% dplyr::select(position, z, p, 
        pop) %>% mutate(region = region, pop = case_when(pop == "bbj-a-52" ~ "EAS", pop == "ieu-b-38" ~ "EUR", pop == "fema" ~ "FEMA"))
    print(temp)
    th <- temp %>% dplyr::group_by(pop) %>% dplyr::arrange(desc(z)) %>% 
        dplyr::slice_head(n = 1) %>% dplyr::ungroup()
    the <- temp %>% filter(pop == "EUR") %>% dplyr::arrange(desc(z)) %>% 
        dplyr::slice_head(n = 1)
    the <- subset(temp, position == the$position)
    thf <- subset(temp, position == subset(th, pop == "FEMA")$position)
    ggplot2::ggplot(temp, ggplot2::aes(x = position, y = -log10(p))) + 
        ggplot2::geom_point(colour="grey") + 
        ggplot2::geom_segment(data = th, 
            aes(x = position, xend = position, y = 0, yend = -log10(p)), 
            colour = "pink") + 
        ggplot2::geom_point(data = th, colour = "red") + 
        ggplot2::geom_segment(data = thf, aes(x = position, xend = position, 
            y = 0, yend = -log10(p)), colour = "grey") + 
        ggplot2::geom_point(data = thf, 
            colour = "blue") + 
        ggplot2::geom_point(data = the, colour = "black") + 
        ggplot2::facet_grid(pop ~ ., scale="free_y") +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
}

6
13
36
43
49
61


i <- 43
pl <- p(names(x$instrument_regions)[i])
pl
ggsave(pl, file=here("results", "fema_example.pdf"), width=4, height=6)
ggsave(pl, file=here("results", "fema_example.png"), width=4, height=6)
