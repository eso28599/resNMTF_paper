#!/bin/bash
#PBS -N increasing_noise
#PBS -m a
#PBS -q medium
#PBS -t 1-100
#PBS -o Results/noise/noise_3v4b/logs/test_job.out
#PBS -e Results/noise/noise_3v4b/logs/test_job.err

export R_LIBS="/home/clustor2/ma/e/eso18/R/x86_64-pc-linux-gnu-library/4.3"
export sim_folder_name=Results/noise/noise_3v4b
export sim=noise
export i=${PBS_ARRAYID}
export I=`echo $i | awk '{printf "%3.3d", $1}'`


cd ${PBS_O_WORKDIR}/${sim_folder_name}/data
#mkdir $I
if [ ! -d "$I" ]; then
  mkdir $I
  cd $I
  for i in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100
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
python3 OtherMethods/methods_p.py ${sim_folder_name} $I ${sim}

#evaluate results in R
Rscript --vanilla eval.r  ${sim_folder_name} $I
