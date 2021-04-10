library(ieugwasr)
library(gwasglue)
library(TwoSampleMR)
library(susieR)
library(ggplot2)
library(dplyr)


susie_overlaps <- function(su1, su2)
{
	l <- list()
	k <- 1
	s1 <- su1$sets$cs
	s2 <- su2$sets$cs
	for(i in 1:length(s1))
	{
		for(j in 1:length(s2))
		{
			if(any(s1[[i]] %in% s2[[j]]))
			{
				ind <- s1[[i]] %in% s2[[j]]
				v <- s1[[i]][ind]
				l[[k]] <- tibble(s1=i, s2=j, v = v)
				k <- k + 1
			}
		}
	}
	l <- bind_rows(l)
	if(nrow(l) > 0)
	{
		l$pip1 <- su1$pip[l$v]
		l$pip2 <- su2$pip[l$v]
		l$piprank1 <- rank(-su1$pip)[l$v] / length(su1$pip)
		l$piprank2 <- rank(-su2$pip)[l$v] / length(su2$pip)
		l <- l %>%
			mutate(piprank=piprank1 + piprank2) %>%
			arrange(piprank)
	}
	return(l)
}


greedy_remove <- function(r, thresh)
{
  diag(r) <- 0
  r <- abs(r)
  flag <- 1
  rem <- c()
  nom <- colnames(r)
  while(flag == 1)
  {
      count <- apply(r, 2, function(x) sum(x >= thresh))
      if(any(count > 0))
      {
          worst <- which.max(count)[1]
          rem <- c(rem, names(worst))
          r <- r[-worst,-worst]
      } else {
          flag <- 0
      }
  }
  return(which(nom %in% rem))
}

independent_instruments <- function(id1, id2, alpha=0.05, method="fdr", ld_thresh=0.05, bfile=NULL, plink=NULL, ld_pop="EUR")
{
	message("Tophits for ", id1)
	p1th <- tophits(id1)
	message("Tophits for ", id2)
	p2th <- tophits(id2)
	message("Looking up reciprocal associations")
	lu <- dplyr::bind_rows(
		associations(p2th$rsid, id1) %>%
			dplyr::select(rsid, chr, position, ea, nea, id, beta, se, p),
		associations(p1th$rsid, id2) %>%
			dplyr::select(rsid, chr, position, ea, nea, id, beta, se, p)
	)
	lu$pvaladjust <- p.adjust(lu$p, method)
	rsid <- unique(lu$rsid)
	message("Found ", length(rsid), " unique variants")
	message("LD pruning")
	ld <- suppressMessages(suppressWarnings(ieugwasr::ld_matrix(rsid, pop=ld_pop, bfile=bfile, plink=plink)))
	ldnom <- rownames(ld) %>% strsplit(., split="_") %>% sapply(., function(x) x[1])
	rem <- greedy_remove(ld, sqrt(ld_thresh))
	to_remove <- ldnom[rem] %>% strsplit(., split="_") %>% sapply(., function(x) x[1])
	to_remove <- c(to_remove, rsid[!rsid %in% ldnom])
	message("Removing ", length(to_remove), " due to LD or absence from reference")
	lu <- subset(lu, !rsid %in% to_remove) %>%
		mutate(sig = pvaladjust < alpha)
	message(nrow(lu), " variants remain")

	# shared_snps <- subset(lu, sig)
	# unique_snps1 <- subset(lu, p > 0.05 & id == id2)
	# unique_snps2 <- subset(lu, p > 0.05 & id == id1)
	return(lu)
}

trans_mr_instruments <- function(id1, id2, radius, pop1=NULL, pop2=NULL, plink=NULL, bfile1=NULL, bfile2=NULL, ld_thresh=0.05, alpha=0.05, mt_method="fdr", prune_bfile=NULL, prune_pop=NULL)
{
	lu <- independent_instruments(id1, id2, alpha, mt_method, ld_thresh, prune_bfile, plink, prune_pop)
	l <- list()
	for(i in 1:nrow(lu))
	{
		message(lu$chr[i], ":", lu$position[i])
		l[[i]] <- finemap_region(chr=lu$chr[i], position=lu$position[i], radius=radius, id1=id1, id2=id2, plink=plink, bfile1=bfile1, bfile2=bfile2, pop1=pop1, pop2=pop2)
	}

	res <- lapply(l, function(x)
	{
		tibble(
			type=ifelse(is.null(x$type), NA, x$type), 
			bestsnp=ifelse(is.null(x$bestsnp), NA, x$bestsnp), 
			cs_overlap=ifelse(is.null(x$cs_overlap), NA, x$cs_overlap)
		)
	}) %>% bind_rows()

	inst <- extract_outcome_data(subset(res, !is.na(bestsnp))$bestsnp, c(id1, id2)) %>% convert_outcome_to_exposure()
	return(list(inst=inst, res=res, info=l))
}

finemap_region <- function(chr, position, radius, id1, id2, plink=NULL, bfile1=NULL, bfile2=NULL, pop1=NULL, pop2=NULL)
{
	r <- paste0(chr, ":", max(0, position-radius), "-", position+radius)
	message("Analysing ", r)
	message("Extracting summary data and LD for ", id1)
	a1 <- ieugwasr::associations(r, id1) %>%
		arrange(position) %>%
		filter(!duplicated(rsid)) %>%
		gwasglue::ieugwasr_to_TwoSampleMR(.)
	a1 <- TwoSampleMR::harmonise_ld_dat(a1, suppressMessages(suppressWarnings(ieugwasr::ld_matrix(a1$SNP, pop=pop1, plink=plink, bfile=bfile1, with_alleles=TRUE))))
	message("Extracting summary data and LD for ", id2)
	a2 <- ieugwasr::associations(r, id2) %>% 
		arrange(position) %>%
		filter(!duplicated(rsid)) %>%
		gwasglue::ieugwasr_to_TwoSampleMR(.)
	a2 <- TwoSampleMR::harmonise_ld_dat(a2, suppressMessages(suppressWarnings(ieugwasr::ld_matrix(a2$SNP, pop=pop2, plink=plink, bfile=bfile2, with_alleles=TRUE))))

	snps <- a1$x$SNP[a1$x$SNP %in% a2$x$SNP]
	message(length(snps), " variants in common")
	index1 <- which(a1$x$SNP %in% snps)
	index2 <- which(a2$x$SNP %in% snps)
	a1$x <- a1$x[index1,]
	a1$ld <- a1$ld[index1, index1]
	a2$x <- a2$x[index2,]
	a2$ld <- a2$ld[index2, index2]

	stopifnot(all(a1$x$SNP == a2$x$SNP))

	message("Running susie for ", id1)
	su1 <- susieR::susie_rss(z=a1$x$beta.exposure / a1$x$se.exposure, R=a1$ld, check_R=FALSE)
	message("Running susie for ", id2)
	su2 <- susieR::susie_rss(z=a2$x$beta.exposure / a2$x$se.exposure, R=a2$ld, check_R=FALSE)

	message("Finding overlaps")
	suo <- susie_overlaps(su1, su2)
	out <- list(
		chr=chr,
		position=position,
		radius=radius,
		a1=a1$x,
		a2=a2$x,
		su1=su1,
		su2=su2,
		suo=suo
	)

	null <- c(is.null(su1$sets$cs), is.null(su2$sets$cs))
	if(null[1] & !null[2])
	{
		out$type <- "pop2"
	} else if(null[2] & !null[1]) {
		out$type <- "pop1"
	} else if(!null[1] & !null[2]) {
		out$type <- "shared"
	} else {
		out$type <- "drop"
	}

	if(out$type == "shared")
	{
		if(nrow(suo) == 0)
		{
			out$cs_overlap <- FALSE
			temp <- inner_join(a1$x, a2$x, by="SNP") %>%
				mutate(pvalrank = rank(pval.exposure.x) / nrow(a1$x) + rank(pval.exposure.y) / nrow(a1$x)) %>%
				arrange(pvalrank)
			out$bestsnp <- temp$SNP[1]
		} else {
			out$cs_overlap <- TRUE
			out$bestsnp <- a1$x$SNP[suo$v[1]]
		}
	}

	if(out$type == "pop1")
	{
		out$cs_overlap <- FALSE
		out$bestsnp <- a1$x %>% arrange(pval.exposure) %>% {.$SNP[1]}
	}
	if(out$type == "pop2")
	{
		out$cs_overlap <- FALSE
		out$bestsnp <- a2$x %>% arrange(pval.exposure) %>% {.$SNP[1]}
	}
	return(out)
}

plot_region <- function(region)
{
	lapply(list(region$a1, region$a2), function(x)
		x %>%
			dplyr::select(pval.exposure, position.exposure, id.exposure)
		) %>% bind_rows() %>%
		ggplot(., aes(x=position.exposure, y=-log10(pval.exposure))) +
		geom_point() +
		facet_grid(id.exposure ~ .)
}




radius=50000
id1 <- "bbj-a-52"
id2 <- "ukb-b-20175"
bfile1="/Users/gh13047/data/ld_files/ldmat/EAS"
bfile2="/Users/gh13047/data/ld_files/ldmat/EUR"
plink="plink"
alpha=0.05
mt_method="fdr"
ld_thresh=0.05
prune_bfile=NULL
prune_pop="EUR"

inst <- trans_mr_instruments(id1, id2, radius, plink=plink, bfile1=bfile1, bfile2=bfile2, ld_thresh=ld_thresh, prune_bfile=bfile1)


d <- associations(inst$res$bestsnp[!is.na(inst$res$bestsnp)], args$ids)
d$padj <- p.adjust(d$p, "fdr")
dplyr::select(d, rsid, id, padj) %>%
	dplyr::group_by(rsid, id) %>%
	dplyr::slice(1) %>%
	pivot_wider(names_from=id, values_from=padj) %>%
	ungroup() %>%
	summarise(s1 = sum(`bbj-a-52` < 0.05 & `ukb-b-20175` > 0.5), s2 = sum(`bbj-a-52` > 0.5 & `ukb-b-20175` < 0.05), sb = sum(`bbj-a-52` < 0.05 & `ukb-b-20175` < 0.05))

inst$inst %>%
	mutate(padj = p.adjust(pval.exposure, "fdr")) %>%
	dplyr::select(rsid=SNP, id=id.exposure, padj) %>%
	dplyr::group_by(rsid, id) %>%
	dplyr::slice(1) %>%
	pivot_wider(names_from=id, values_from=padj) %>%
	ungroup() %>%
	summarise(s1 = sum(`bbj-a-52` < 0.05 & `ukb-b-20175` > 0.5), s2 = sum(`bbj-a-52` > 0.5 & `ukb-b-20175` < 0.05), sb = sum(`bbj-a-52` < 0.05 & `ukb-b-20175` < 0.05))

table(d$rsid %in% inst$inst$SNP)

d3 <- independent_instruments(args$ids[1], args$ids[2])

d4 <- associations(d3$rsid, args$ids)
d4$padj <- p.adjust(d4$p, "fdr")
dplyr::select(d4, rsid, id, padj) %>%
	dplyr::group_by(rsid, id) %>%
	dplyr::slice(1) %>%
	pivot_wider(names_from=id, values_from=padj) %>%
	ungroup() %>%
	summarise(s1 = sum(`bbj-a-52` < 0.05 & `ukb-b-20175` > 0.5), s2 = sum(`bbj-a-52` > 0.5 & `ukb-b-20175` < 0.05), sb = sum(`bbj-a-52` < 0.05 & `ukb-b-20175` < 0.05))



inst$res

which(inst$res$cs_overlap)


plot_region(inst$info[[100]])
susie_plot(inst$info[[100]]$su1, "PIP")
dev.new()
susie_plot(inst$info[[100]]$su2, "PIP")

table(inst$res$cs_overlap)

dim(inst$res)
table(inst$res$type)

saveRDS(inst, file="sbp_eas_eur.rds")

plot_region(l[[200]])
plot_region(l[[207]])

susie_plot(l[[2]]$su1, "PIP")
susie_plot(l[[2]]$su2, "PIP")

