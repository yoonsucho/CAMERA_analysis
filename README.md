# Trans ancestral design for Mendelian Randomization

Develop trans-ancestral model for MR

## Background
Emerging availability of large scale GWAS in non-European populations has enabled `trans-ancestral' studies that are conducted using more than one population. Trans-ancestral study design have been suggested to overcome existing issues with limited sample sizes available for many diseases. 

There have been a few studies that introduced trans-ancestral design to MR. However, no guideline exist for conducting and reporting trans-ancestral MR so far. Here, we simulated possible scenarios for trans-ancestral MR to provide a guideline for the trans-ancestral MR study design.

---

## To do:

1. Simulations for population differences
2. Pleiotropy model
3. Simulations for pleiotropy
4. Making an R package
5. Apply to real data


### Population simulation questions

Power and false discovery rates to detect differences in causal effects in two populations by

- Different sample sizes
- Different allele frequencies
- Different correlations between X1 and X2
- Different heritabilities of X1 and X2
- Different values of biv1 and biv2
- Different thresholds of instrument p-values (introduces winner's curse)

---

## Pleiotropy analysis outline

### Instrument generation

Exposure and outcome data is available in two populations. e.g. SBP and CHD

1. Tophits from pop1
2. Tophits from pop2
3. Lookup tophits from pop1 in pop2 - X1 = anything in pop1 that has p > 0.05 in pop2
4. Lookup tophits from pop2 in pop1 - X2 = anything in pop2 that has p > 0.05 in pop1
5. Create Xb. Anything from pop1 or pop2 that has p < 0.05 / nsnp and the effect is in the same direction in pop2 or pop1 is to be kept
6. Finemap SNPs in pop1 - get credible sets in each region
7. Finemap SNPs in pop2 - get credible sets in each region
8. Find overlapping credible sets between the populations - define the best SNP by ranking the shared SNPs by p-value in each pop, summing the rank and keeping the lowest value
9. If there are multiple credible sets in the region, just keep one shared effect and discard everything else
10. Meta-analyse the effects (IVW) to get Xb
11. Create a data from of exposure data - X1, Xb, X2

### Outcome effects

1. Lookup all SNPs in formatted exposure data in the outcome datasets for each population
2. Update instrument dataframe to add Y1 and Y2 columns

### Analysis

Use SEM models

## Simulation study

Check if the pleiotropy estimation method works under different scenarios

- Directional pleiotropy
- Pleiotropy correlated with SNP-exposure effect
- False negatives in the SNP1-exposure2 estimates

