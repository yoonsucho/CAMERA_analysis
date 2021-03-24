# Trans ancestral design for Mendelian Randomization
Develop trans-ancestral model for MR

## Background
Emerging availability of large scale GWAS in non-European populations has enabled `trans-ancestral' studies that are conducted using more than one population. Trans-ancestral study design have been suggested to overcome existing issues with limited sample sizes available for many diseases. 

There have been a few studies that introduced trans-ancestral design to MR. However, no guideline exist for conducting and reporting trans-ancestral MR so far. Here, we simulated possible scenarios for trans-ancestral MR to provide a guideline for the trans-ancestral MR study design.

### Can we perform two sample MR using the samples from different ancestry?
MR studies require ethnically homogenous samples. The following scenario can be developed:
- Size and direction of genetic effect can be different across the ancestries
- Causal structures are different across the ancestries.
- Minor allele frequencies are different
- Ancestry specific SNPs

### Simulation
1. MR of exposure in ASN against outcome in EUR or vice versa
-	Set different MAF, Effect size (e.g. How big the difference should be?)
- Suggested model: Y1,Y2 ~ X1 + X2, where X1 and Y1 is the exposure and the outcome in ASN, and X2 and Y2 is the expousre and the outcome in EUR, respectively.
2. MR of exposure and outcome in ASN using EUR discovery SNPs
3. Using SEM model to get an estimate of pleiotropy 
4. Using fine-mapping to help improve analyses


