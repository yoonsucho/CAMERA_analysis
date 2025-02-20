# Whole genome LD loss

Motivation for a cross-ancestry MR method needs to be expanded. When you have GWAS discovery in Europeans then it's tempting to estimate the MR effect of an exposure on a non-European outcome.

$$
\beta_{iv} = \beta_{gy} / \beta_{gx}
$$

If the causal variant for the exposure isn't identified then there is a risk that in the non-European population $\beta_{gx}$ is underestimated due to incomplete LD.

$$
\begin{aligned}
\beta_{iv,cross} & = \beta_{gy,non-eur} / \beta_{gx,eur} \\
& = \beta_{iv} \beta_{gx,non-eur} / \beta_{gx,eur} \\
& = \beta_{iv} \beta_{gx,eur} r^2_{eur,non-eur} / \beta_{gx,eur} \\
& = \beta_{iv} r^2_{eur,non-eur}
\end{aligned}
$$

Alternatively, if the effect estimate is based on re-estimated $\beta_{gx}$ and $\beta_{gy}$ in the non-European population then even if there is incomplete LD that should cancel out. However, the larger the loss in LD the weaker the non-Eur IV will be, and similarly if there is ascertainment for more intermediate allele frequencies in Europeans then weaker IV in non-Eur.

## Questions

1. What is the distribution of LD-loss for widely used European derived instruments?
    - Degree of bias in non-Eur
    - Loss of power in non-Eur
    - Are Eur instruments ascertained for differential LD regions?
2. What is the distribution of AF reduction for widely used European derived instruments?
    - Loss of power in non-Eur
    - Are Eur instruments ascertained for differential LD regions?

## Methods

Use 1000 genomes reference panel, LD and allele frequency info for 5 super populations.

1. Select a causal variant from European sample. Give it a p-value
2. For all SNPs in the region, what is the probability that it is selected as the top hit?
3. For each SNP that has non-zero probability of being selected, what is the LD loss?
4. Integrating over 2 and 3, what is the LD loss distribution for each SNP?

## Other simulations

Need to organise the simulations that were performed.

1. Using maxz to improve instrument selection
2. Using SEM to jointly estimate effects

## Coding notes

- Using snakemake to pipeline the analysis
- All analysis in R so using renv to control environment
- Prototyping locally but distributing parallel code to bluecrystal4

To run

Setup `config.json` like:

```json
{
    "ld_data_dir": "somewhere",
    "outdir": "somewhere"
}
```

1. Install `snakemake` e.g. 
    ```
    conda env create -f env.yaml
    conda activate snakemakemin
    ```
2. Install renv package and then restore using `renv::restore()`


```
snakemake
```

On bc4

```
snakemake -prk \
-j 100 \
--cluster-config bc4-cluster.json \
--cluster "sbatch \
  --job-name={cluster.name} \
  --partition={cluster.partition} \
  --nodes={cluster.nodes} \
  --ntasks-per-node={cluster.ntask} \
  --cpus-per-task={cluster.ncpu} \
  --time={cluster.time} \
  --mem={cluster.mem} \
  --account={cluster.account} \
  --output={cluster.output}"
```


