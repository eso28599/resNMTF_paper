#!/bin/bash
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -l walltime=05:00:00
#PBS -N other_results
#PBS -o /rds/general/user/eso18/home/resNMTF_paper/RealData/logs/run_others.out
#PBS -e /rds/general/user/eso18/home/resNMTF_paper/RealData/logs/run_others.err


cd resNMTF_paper
# python methods
# eval "$(~/miniforge3/bin/conda shell.bash hook)" 
# conda activate resnmtf_env
# python3 RealData/3sources/3sources_python_analysis.py
# python3 RealData/bbcsport/bbc_python_analysis.py
# # python3 RealData/single_cell/single_cell_python_analysis.py
# # python3 RealData/3cancers/3cancers_python_analysis.py
# conda deactivate


# remove bad package installation
# rm -rf /rds/general/user/eso18/home/R/x86_64-pc-linux-gnu-library/4.4/00LOCK-ggplot2
# rm -rf /rds/general/user/eso18/home/R/x86_64-pc-linux-gnu-library/4.4/00LOCK-biclust
# rm -rf /rds/general/user/eso18/home/R/x86_64-pc-linux-gnu-library/4.4/biclust
# rm -rf /rds/general/user/eso18/home/R/x86_64-pc-linux-gnu-library/4.4/ggplot2

# R based other methods
module purge
ml tools/prod #always have to be loaded
module load  R/4.4.2-gfbf-2024a 
# Rscript -e 'install.packages("biclust", type="source", repos="https://cran.ma.imperial.ac.uk/")'

export R_LIBS_USER=/rds/general/user/eso18/home/R/x86_64-pc-linux-gnu-library/4.4 
Rscript --vanilla RealData/3sources/3sources_other_results.r

# combine results

# Rscript --vanilla RealData/combine_dis_metrics.r 3sources 3