#!/bin/bash
#PBS -l select=1:ncpus=5:mem=30gb
#PBS -l walltime=08:00:00
#PBS -N 3sources_stab
#PBS -o /rds/general/user/eso18/home/resNMTF_paper/RealData/bbcsport/logs
#PBS -e /rds/general/user/eso18/home/resNMTF_paper/RealData/bbcsport/logs


module load  R/4.4.2-gfbf-2024a 
export R_LIBS_USER=/rds/general/user/eso18/home/R/x86_64-pc-linux-gnu-library/4.4  
cd resNMTF_paper
Rscript --vanilla RealData/bbcsport/bbc_stab.r 1
