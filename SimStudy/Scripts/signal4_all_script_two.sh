#!/bin/bash
#PBS -N issvd_4all_res4
#PBS -m a
#PBS -q medium
#PBS -o ../Results/signal/signal4_all/logs/results.out
#PBS -e ../Results/signal/signal4_all/logs/results.err

export R_LIBS="/home/clustor4/ma/e/eso18/R/x86_64-pc-linux-gnu-library/4.4"
export sim_folder_name=signal/signal4_all
export sim=signal


cd ${PBS_O_WORKDIR}/../Results/${sim_folder_name}
ls data/ > repeat_folders.txt
#move back into original folder
cd ${PBS_O_WORKDIR}/..


#compile results
Rscript -e "rmarkdown::render('results.rmd')" Results/${sim_folder_name}  Results/${sim_folder_name}/repeat_folders.txt

mv results.pdf Results/${sim_folder_name}
