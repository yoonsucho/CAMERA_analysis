library(ieugwasr)
library(gwasglue)
library(TwoSampleMR)
library(susieR)
library(ggplot2)
library(dplyr)

MultiAncestrySummarySet <- R6::R6Class("MultiAncestrySummarySet", list(
  output = list(),
  exposure_ids=NULL,
  outcome_ids=NULL,
  pops=NULL,
  bfiles=NULL,
  plink=NULL,
  raw_instruments=NULL,
  instrument_regions = NULL,
  ld_matrices = NULL,
  susie_results = NULL,
  paintor_results = NULL,
  expected_replications = NULL,



  # Methods

  initialize = function(exposure_ids, outcome_ids, pops=NULL, bfiles=NULL, plink=NULL)
  {
  	self$exposure_ids <- exposure_ids
  	self$outcome_ids <- outcome_ids
  	self$pops <- pops
  	self$bfiles <- bfiles
  	self$plink <- plink
  },

  # e.g. could just use mv_extract_exposures
  # - for each exposure_id it identifies the instruments
  # - then it tries to make these independent
  # - then looks up the same set of instruments in all exposures
  extract_instruments = function(exposure_ids=self$exposure_ids, ...) {
  	self$raw_instruments <- mv_extract_exposures(exposure_ids, ...)
  },

  # Here the idea is that pop1 and pop2 might share an instrument, but the tophit for pop1 is not the causal variant
  # Hence, in pop1 it is in LD with the causal variant but not in pop2
  # So we extract a region around each instrument (e.g. 50kb)
  # Search for a SNP that is best associated in both pop1 and pop2
  extract_instrument_regions = function(raw_instruments=self$raw_instruments, bfiles=self$bfiles, )
  {
	
	ldflag <- !(is.null(bfiles) & is.null(pops))

	if(!is.null(bfiles))
	{
		stopifnot(length(bfiles) == length(ids))
	}
	if(!is.null(pops))
	{
		stopifnot(length(pops) == length(ids))
	}
	r <- paste0(chr, ":", max(0, position-radius), "-", position+radius)
	message("Analysing ", r)

	a <- lapply(1:length(ids), function(i)
	{
		message("Extracting summary data and LD for ", ids[i])
		a1 <- ieugwasr::associations(r, ids[i]) %>%
			dplyr::arrange(position) %>%
			dplyr::filter(!duplicated(rsid)) %>%
			gwasglue::ieugwasr_to_TwoSampleMR(.)
		if(ldflag)
		{
			a1 <- TwoSampleMR::harmonise_ld_dat(a1, suppressMessages(suppressWarnings(ieugwasr::ld_matrix(a1$SNP, pop=pops[i], plink=plink, bfile=bfiles[i], with_alleles=TRUE))))
			return(a1)
		} else {
			return(list(x=a1))
		}
	})

	snps <- Reduce(intersect, lapply(a, function(x) x$x$SNP))
	a <- lapply(a, function(x)
	{
		index1 <- which(x$x$SNP %in% snps)
		x$x <- x$x[index1,]
		if(ldflag)
		{
			x$ld <- x$ld[index1, index1]
		}
		return(x)
	})
	if(length(ids) > 1)
	{
		stopifnot(all(a[[1]]$x$SNP == a[[2]]$x$SNP))
	}
	return(a)  	
  },

  # to build on extract_instrument_regions we can do finemapping
  # If we want to do fine mapping we need to get an LD matrix for the whole region (for each population)
  # We then need to harmonise the LD matrix to the summary data, and the summary datasets to each other
  regional_ld_matrices = function(),

  # for each instrument region + ld matrix we can perform susie finemapping
  # do this independently in each population
  # find credible sets that overlap - to try to determine the best SNP in a region to be used as instrument
  susie_finemap_regions = function(),

  # PAINTOR allows finemapping jointly across multiple populations
  # returns a posterior probability of inclusion for each SNP
  # Choose the variant with highest posterior probability and associations in each exposure
  paintor_finemap_regions = function(),

  # Once a set of instruments is chosen we can ask
  # - what fraction of those primarily identified in pop1 replicate in pop2?
  # - what fraction of those primarily identified in pop1 have the same sign as in pop2?
  # - compare these to what is expected by chance under the hypothesis that the effect estimates are the same
  # this will provide some evidence for whether lack of replication is due to (e.g.) GxE
  proportion_instrument_overlap = function(b_disc, b_rep, se_disc, se_rep, alpha)
  {
	p_sign <- pnorm(-abs(b_disc) / se_disc) * pnorm(-abs(b_disc) / se_rep) + ((1 - pnorm(-abs(b_disc) / se_disc)) * (1 - pnorm(-abs(b_disc) / se_rep)))
	p_sig <- pnorm(-abs(b_disc) / se_disc + qnorm(alpha / 2)) + (1 - pnorm(-abs(b_disc) / se_disc - qnorm(alpha / 2)))
	p_rep <- pnorm(abs(b_rep)/se_rep, lower.tail=FALSE)
	message("Total: ", length(b_disc))
	message("Expected sign: ", sum(p_sign, na.rm=TRUE))
	message("Observed sign: ", sum(sign(b_disc) == sign(b_rep)))
	message("Expected sig: ", sum(p_sig, na.rm=TRUE))
	message("Observed sig: ", sum(p_rep < alpha, na.rm=TRUE))
	return(list(sum(p_sign, na.rm=TRUE), sum(p_sig, na.rm=TRUE), data_frame(sig=p_sig, sign=p_sign)))
  },

  # Once a set of instruments is defined then extract the outcome data
  # Could use 
  # - raw results from extract_instruments
  # - maximised associations from extract_instrument_regions
  # - finemapped hits from susie_finemap_regions
  # - finemapped hits from paintor_finemap_regions
  extract_outcome_data = function(),

  # Perform basic SEM analysis of the joint estimates in multiple ancestries
  perform_basic_sem = function()


  # Plots
  # Plot regional associations for each instrument and each population
  # Plot finemapping results against regional association data

))









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

#' Analyse regional data with PAINTOR
#'
#' @param regiondata Output from extract_regional_data
#' @param PAINTOR Path to PAINTOR executable. Default="PAINTOR" 
#' @param workdir Working directory. Default=tempdir()
#'
#' @export
#' @return Results table with posterior inclusion probabilities
run_PAINTOR <- function(regiondata, PAINTOR="PAINTOR", workdir=tempdir())
{
	nid <- length(regiondata)
	snps <- regiondata[[1]]$x$SNP
	zs <- lapply(1:nid, function(i)
	{
		o <- list()
		o[[paste0("ZSCORE.P", i)]] <- regiondata[[i]]$x$beta.exposure / regiondata[[i]]$x$se.exposure
		return(as_tibble(o))
	}) %>% bind_cols()

	locus <- tibble(CHR = regiondata[[1]]$x$chr.exposure, POS = regiondata[[1]]$x$position.exposure, RSID = regiondata[[1]]$x$SNP) %>%
		bind_cols(., zs)
	write.table(locus, file=file.path(workdir, "Locus1"), row=F, col=T, qu=F)
	lapply(1:nid, function(i)
	{
		write.table(regiondata[[i]]$ld, file=file.path(workdir, paste0("Locus1.LD", i)), row=FALSE, col=FALSE, qu=FALSE)
	})
	write.table(tibble(null=rep(1, length(snps))), file=file.path(workdir, "Locus1.annotations"), row=F, col=T, qu=F)

	write.table(c("Locus1"), file=file.path(workdir, "input.files"), row=F, col=F, qu=F)

	LDname <- paste(paste0("LD", 1:nid), collapse=",")
	Zhead <- paste(names(zs), collapse=',')
	glue("{PAINTOR} -input {file.path(workdir,'input.files')} -Zhead {Zhead} -LDname {LDname} -in {workdir} -out {workdir} -mcmc -annotations null") %>%
		system()
	res <- data.table::fread(file.path(workdir, "Locus1.results"))
	unlink(workdir)
	return(res)
}


#' Plot PAINTOR results
#'
#' @param res Output from run_PAINTOR
#'
#' @export
#' @return plot
plot_PAINTOR <- function(res)
{
	res %>% pivot_longer(cols=c(ZSCORE.P1, ZSCORE.P2, Posterior_Prob)) %>%
		ggplot(., aes(x=POS, y=value)) +
		geom_point(aes(colour=name)) +
		facet_grid(name ~ ., scale="free_y")
}

#' Analyse regional data with MsCAVIAR
#'
#' @param regiondata Output from extract_regional_data
#' @param MsCAVIAR Path to MsCAVIAR executable. Default="MsCAVIAR" 
#' @param workdir Working directory. Default=tempdir()
#'
#' @export
#' @return Results table with posterior inclusion probabilities
run_MsCAVIAR <- function(regiondata, n, MsCAVIAR="MsCAVIAR", workdir=tempdir())
{
	nid <- length(regiondata)
	snps <- regiondata[[1]]$x$SNP
	lapply(1:nid, function(i)
	{
		regiondata[[i]]$x %>% 
			dplyr::mutate(z=beta.exposure/se.exposure) %>%
			dplyr::select(SNP, z) %>%
			write.table(., file=file.path(workdir, paste0("z.", i)), row=F, col=F, qu=F)
	})

	lapply(1:nid, function(i)
	{
		write.table(regiondata[[i]]$ld, file=file.path(workdir, paste0("ld.", i)), row=FALSE, col=FALSE, qu=FALSE)
	})

	thisdir <- getwd()
	setwd(workdir)
	write.table(paste0("ld.", 1:nid), file="ldfiles.txt", row=F, col=F, qu=F)
	write.table(paste0("z.", 1:nid), file="zfiles.txt", row=F, col=F, qu=F)
	glue("{MsCAVIAR} -l ldfiles.txt -z zfiles.txt -n {paste(n, collapse=',')} -o results") %>% system()
	res <- data.table::fread("results_post.txt")
	zs <- lapply(1:nid, function(i)
	{
		out <- list()
		out[[paste0("ZSCORE.P", i)]] <- regiondata[[i]]$x$beta.exposure / regiondata[[i]]$x$se.exposure
		as_tibble(out)
	}) %>% bind_cols()
	stopifnot(all(res$SNP_ID == regiondata[[1]]$x$SNP))
	names(res)[1] <- c("SNP")
	res <- bind_cols(res, zs) %>%
		mutate(CHR = regiondata[[1]]$x$chr.exposure, POS=regiondata[[1]]$x$position.exposure)
	setwd(thisdir)
	return(res)
}

#' Plot MsCAVIAR results
#'
#' @param res Output from run_MsCAVIAR
#'
#' @export
#' @return plot
plot_MsCAVIAR <- function(res)
{
	res %>% pivot_longer(cols=c(ZSCORE.P1, ZSCORE.P2, Prob_in_pCausalSet)) %>%
		ggplot(., aes(x=POS, y=value)) +
		geom_point(aes(colour=name)) +
		facet_grid(name ~ ., scale="free_y")
}

