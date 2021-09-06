
```
module add languages/anaconda3/5.2.0-tflow-1.11
snakemake -prk \
-j 100 \
--cluster-config bc4-cluster.json \
--cluster "sbatch \
  --job-name={cluster.name} \
  --partition={cluster.partition} \
  --ntasks-per-node={cluster.ntask} \
  --cpus-per-task={cluster.ncpu} \
  --time={cluster.time} \
  --mem={cluster.mem} \
  --output={cluster.output}"
```


