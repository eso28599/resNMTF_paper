#!/bin/bash
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -l walltime=01:00:00
#PBS -N 3sources
#PBS -o /rds/general/user/eso18/home/resNMTF_paper/RealData/3sources/logs/run_others.out
#PBS -e /rds/general/user/eso18/home/resNMTF_paper/RealData/3sources/logs/run_others.err


module load  R/4.4.2-gfbf-2024a 
export R_LIBS_USER=/rds/general/user/eso18/home/R/x86_64-pc-linux-gnu-library/4.3  
cd resNMTF_paper

# # # python scripts
# eval "$(~/miniforge3/bin/conda shell.bash hook)" to load conda
# conda activate resnmtf_env
# # python3 RealData/3sources/3sources_python_analysis.py
# python3 RealData/single_cell/single_cell_python_analysis.py

Rscript --vanilla RealData/3sources/3sources_other_results.r

