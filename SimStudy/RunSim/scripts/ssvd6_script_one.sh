#!/bin/bash
#PBS -N issvd_param
#PBS -m a
#PBS -q medium
#PBS -t 1-100
#PBS -o Results/issvd_param6/logs/test_job.out
#PBS -e Results/issvd_param6/logs/test_job.err

export R_LIBS="/home/clustor4/ma/e/eso18/R/x86_64-pc-linux-gnu-library/4.4"
# export R_LIBS="/usr/local/lib/R/site-library"
export sim_folder_name=issvd_param6
export sim=issvd
export i=${PBS_ARRAYID}
export I=`echo $i | awk '{printf "%3.3d", $1}'`


cd ${PBS_O_WORKDIR}/../Results/${sim_folder_name}/data
#mkdir $I
if [ ! -d "$I" ]; then
  mkdir $I
  cd $I
  for i in {0..21}
  do
    mkdir issvd_$i
  done
fi

#move back into original folder
cd ${PBS_O_WORKDIR}/..

#generate data
Rscript --vanilla data_gen.r  Results/${sim_folder_name} $I

#analyse in python
python3 OtherMethods/issvd_param_p.py Results/${sim_folder_name} $I

#evaluate results in R
Rscript --vanilla eval_ssvd_param.r  Results/${sim_folder_name} $I
