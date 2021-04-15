rm(list=ls())
source("trans_instruments.r")

# Download LD reference data from:
# http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz

# Need to specify path to plink executable also
# These exposures are for CRP
# Outcomes are CHD (but should consider choosing other outcomes for EAS and AFR, or find other traits)
x <- MultiAncestrySummarySet$new(
	exposure_ids=c("ukb-e-30710_EAS", "ukb-e-30710_AFR"), 
	outcome_ids=c("ieu-a-7", "bbj-a-109"), 
	bfiles=c("/Users/gh13047/data/ld_files/ldmat/EAS", "/Users/gh13047/data/ld_files/ldmat/AFR"), 
	plink="plink", 
	radius=50000, 
	clump_pop="AFR"
)

# Extract instruments for each exposure, harmonise etc
x$extract_instruments()
x$instrument_raw %>% as.data.frame

# Estimate what fraction of the raw instruments we'd expect to replicate
x$estimate_instrument_specificity(x$instrument_raw)

# Get regions around each instrument
x$extract_instrument_regions()

# Find region-selected instruments by searching across region for best associations
x$scan_regional_instruments()

# Plot regions and their region-selected instruments
x$plot_regional_instruments(instruments=x$instrument_maxz)

# Estimate what fraction of the region-selected instruments we'd expect to replicate
x$estimate_instrument_specificity(x$instrument_maxz)

# obtain LD matrices
x$regional_ld_matrices()

# TODO: Run SuSiE fine mapping
# This should create a new version of x$instrument_maxz and x$instrument_raw (e.g. called x$instrument_susie) which selects the best SNPs based on finemapping
x$susie_fine_mapping()

# TODO: Run PAINTOR fine mapping
# This should create a new version of x$instrument_maxz and x$instrument_raw (e.g. called x$instrument_paintor) which selects the best SNPs based on finemapping
x$paintor_fine_mapping()

# TODO: get outcome data for a particular set of instruments
x$extract_outcome_data()

# TODO: run SEM model
x$perform_basic_sem()












y <- x
x$import(y)

x <- MultiAncestrySummarySet$new()
x$import(y)

