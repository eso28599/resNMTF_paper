#!/bin/bash
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -l walltime=03:00:00
#PBS -N bbc_other
#PBS -o /rds/general/user/eso18/home/resNMTF_paper/RealData/logs/run_bbc.out
#PBS -e /rds/general/user/eso18/home/resNMTF_paper/RealData/logs/run_bbc.err


cd resNMTF_paper
# python methods
eval "$(~/miniforge3/bin/conda shell.bash hook)" 
conda activate resnmtf_env
# python3 RealData/bbcsport/bbc_python_analysis.py
conda deactivate

# R based other methods
module purge
ml tools/prod #always have to be loaded
module load  R/4.4.2-gfbf-2024a 

export R_LIBS_USER=/rds/general/user/eso18/home/R/x86_64-pc-linux-gnu-library/4.4 
Rscript --vanilla RealData/other_results.r bbcsport

# combine results

# Rscript --vanilla RealData/combine_dis_metrics.r 3sources 3