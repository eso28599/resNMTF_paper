#!/bin/bash
#PBS -N issvd_param_res
#PBS -m a
#PBS -q medium
#PBS -o ../Results/issvd_param6/logs/results.out
#PBS -e ../Results/issvd_param6/logs/results.err

export R_LIBS="/home/clustor4/ma/e/eso18/R/x86_64-pc-linux-gnu-library/4.4"
export sim_folder_name=issvd_param6
export sim=issvd


cd ${PBS_O_WORKDIR}/../Results/${sim_folder_name}
ls data/ > repeat_folders.txt
#move back into original folder
cd ${PBS_O_WORKDIR}/..


#compile results
Rscript -e "rmarkdown::render('issvd_results.rmd')" Results/${sim_folder_name}  Results/${sim_folder_name}/repeat_folders.txt

mv results.pdf ${sim_folder_name}
