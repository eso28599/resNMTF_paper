#!/bin/bash
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -l walltime=01:00:00
#PBS -N combine_res
#PBS -o /rds/general/user/eso18/home/resNMTF_paper/RealData/logs/run_combine.out
#PBS -e /rds/general/user/eso18/home/resNMTF_paper/RealData/logs/run_combine.err

module purge
ml tools/prod #always have to be loaded
module load  R/4.4.2-gfbf-2024a 


export R_LIBS_USER=/rds/general/user/eso18/home/R/x86_64-pc-linux-gnu-library/4.3  
cd resNMTF_paper

Rscript --vanilla RealData/combine_dis_metrics.r 3sources 3

