#!/bin/bash

#SBATCH --job-name=ld_loss
#SBATCH --partition=mrcieu,cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --mem=10G
#SBATCH --account=smed001801
#SBATCH --array=1-500

if [[ -z "${SLURM_ARRAY_TASK_ID}" ]]
then
    i=$1
else
    i=${SLURM_ARRAY_TASK_ID}
fi

region=`head -n $i data/region.txt | tail -n 1 | cut -f 4 -d " "`
srun Rscript whole_genome_ld_loss.r ${region}
