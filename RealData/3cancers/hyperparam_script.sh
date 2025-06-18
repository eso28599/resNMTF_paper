#!/bin/bash
#PBS -l select=1:ncpus=5:mem=10gb
#PBS -l walltime=24:00:00
#PBS -N 3cancers_hyperparam_parallel
#PBS -J 1-123
#PBS -o /rds/general/user/eso18/home/resNMTF_paper/RealData/3cancers/logs
#PBS -e /rds/general/user/eso18/home/resNMTF_paper/RealData/3cancers/logs

module purge
ml tools/prod #always have to be loaded
module load  R/4.4.2-gfbf-2024a 
export R_LIBS_USER=/rds/general/user/eso18/home/R/x86_64-pc-linux-gnu-library/4.3  

cd resNMTF_paper
Rscript --vanilla RealData/3cancers/hyperparam.r $PBS_ARRAY_INDEX
