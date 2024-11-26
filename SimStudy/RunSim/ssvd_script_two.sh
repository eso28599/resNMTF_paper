#!/bin/bash
#PBS -N issvd_param_res
#PBS -m a
#PBS -q medium
#PBS -o issvd_param/logs/results.out
#PBS -e issvd_param/logs/results.err

export R_LIBS="/home/clustor4/ma/e/eso18/R/x86_64-pc-linux-gnu-library/4.4"
export sim_folder_name=issvd_param
export sim=issvd


cd ${PBS_O_WORKDIR}/${sim_folder_name}
ls data/ > repeat_folders.txt
#move back into original folder
cd ${PBS_O_WORKDIR}


#compile results
Rscript -e "rmarkdown::render('results.rmd')" ${sim_folder_name}  ${sim_folder_name}/repeat_folders.txt

mv results.pdf ${sim_folder_name}
