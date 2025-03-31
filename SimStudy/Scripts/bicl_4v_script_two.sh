#!/bin/bash
#PBS -N inc_bicl_res_4v
#PBS -m a
#PBS -q medium
#PBS -o ../Results/bicl/bicl_4v/logs/results.out
#PBS -e ../Results/bicl/bicl_4v/logs/results.err

export R_LIBS="/home/clustor2/ma/e/eso18/R/x86_64-pc-linux-gnu-library/4.3"
export sim_folder_name=Results/bicl/bicl_4v
export sim=bicl

cd ${PBS_O_WORKDIR}/../${sim_folder_name}
ls data/ > repeat_folders.txt
#move back into original folder
cd ${PBS_O_WORKDIR}/..


#compile results
Rscript -e "rmarkdown::render('results.rmd')" ${sim_folder_name}  ${sim_folder_name}/repeat_folders.txt

mv results.pdf ${sim_folder_name}
