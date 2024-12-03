#!/bin/bash
#PBS -N increasing_indiv
#PBS -m a
#PBS -q medium
#PBS -t 1-100
#PBS -o Results/indiv_3v4b/logs/test_job.out
#PBS -e Results/indiv_3v4b/logs/test_job.err

export R_LIBS="/home/clustor2/ma/e/eso18/R/x86_64-pc-linux-gnu-library/4.3"
export sim_folder_name=Results/indiv_3v4b
export sim=indiv
export i=${PBS_ARRAYID}
export I=`echo $i | awk '{printf "%3.3d", $1}'`


cd ${PBS_O_WORKDIR}/${sim_folder_name}/data
#mkdir $I
if [ ! -d "$I" ]; then
  mkdir $I
  cd $I
  for i in 50 200 300 500 1000
  do
    mkdir res_nmtf_$i
    mkdir gfa_$i
    mkdir issvd_$i
    mkdir nmtf_$i
  done
fi

#move back into original folder
cd ${PBS_O_WORKDIR}

#generate data
Rscript --vanilla data_gen.r  ${sim_folder_name} $I

#now analyse in R
Rscript --vanilla methods_r.r  ${sim_folder_name} $I

#analyse in python
python3 OtherMethods/indiv_methods_p.py ${sim_folder_name} $I

#evaluate results in R
Rscript --vanilla eval.r  ${sim_folder_name} $I
