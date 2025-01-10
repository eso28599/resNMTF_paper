#!/bin/bash
#PBS -N inc_indiv_3v4b_res
#PBS -m a
#PBS -q medium
#PBS -o ../Results/indiv/indiv_3v4b/logs/results.out
#PBS -e ../Results/indiv/indiv_3v4b/logs/results.err

export R_LIBS="/home/clustor2/ma/e/eso18/R/x86_64-pc-linux-gnu-library/4.3"
export sim_folder_name=Results/indiv_3v4b
export sim=indiv

cd ${PBS_O_WORKDIR}/../${sim_folder_name}
ls data/ > repeat_folders.txt
#move back into original folder
cd ${PBS_O_WORKDIR}/..


#compile results
Rscript -e "rmarkdown::render('results.rmd')" ${sim_folder_name}  ${sim_folder_name}/repeat_folders.txt

mv results.pdf ${sim_folder_name}
