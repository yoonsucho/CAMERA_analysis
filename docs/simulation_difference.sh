#!/bin/bash
#PBS -l select=1:ncpus=16:mem=100000M
#PBS -l walltime=12:00:00
#PBS -N sim_diff


module add lang/r/4.0.2-gcc
Rscript /work/yc16575/mr.trans/script/simulation_difference.R

