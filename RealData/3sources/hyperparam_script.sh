#!/bin/bash
#PBS -l select=1:ncpus=5:mem=10gb
#PBS -l walltime=02:00:00
#PBS -N 3sources_hyperparam_parallel
#PBS -J 1-123
#PBS -o /rds/general/user/eso18/home/resNMTF_paper/RealData/3sources/logs
#PBS -e /rds/general/user/eso18/home/resNMTF_paper/RealData/3sources/logs

module purge
ml tools/prod #always have to be loaded
module load  R/4.4.2-gfbf-2024a 
export R_LIBS_USER=/rds/general/user/eso18/home/R/x86_64-pc-linux-gnu-library/4.4  
eval "$(~/miniforge3/bin/conda shell.bash hook)" 
conda activate resnmtf_env

cd resNMTF_paper
Rscript --vanilla RealData/3sources/3sources_hyperparam.r $PBS_ARRAY_INDEX
