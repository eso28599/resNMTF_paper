#!/bin/bash
#PBS -N issvd_data_gen
#PBS -m a
#PBS -q medium
#PBS -t 1-69
#PBS -o ../Results/issvd_data_gen/case2/logs/test_job.out
#PBS -e ../Results/issvd_data_gen/case2/logs/test_job.err

export R_LIBS="/home/clustor4/ma/e/eso18/R/x86_64-pc-linux-gnu-library/4.4"
export sim_folder_name=Results/issvd_data_gen/case2
export sim=issvd_data
export i=${PBS_ARRAYID}
export I=`echo $i | awk '{printf "%3.3d", $1}'`


cd ${PBS_O_WORKDIR}/../${sim_folder_name}/data
#mkdir $I
if [ ! -d "$I" ]; then
  mkdir $I
  cd $I
  for i in {1..10}
  do
    mkdir res_nmtf_$i
    mkdir gfa_$i
    mkdir issvd_$i
    mkdir nmtf_$i
  done
fi

#move back into original folder
cd ${PBS_O_WORKDIR}/..

#generate data
Rscript --vanilla issvd_data_gen.r  ${sim_folder_name} $I "two"

# # analyse in r
# Rscript --vanilla methods_r.r  ${sim_folder_name} $I

#analyse in python
python3 OtherMethods/methods_p.py ${sim_folder_name} $I ${sim}

#evaluate results in R
Rscript --vanilla eval.r  ${sim_folder_name} $I
