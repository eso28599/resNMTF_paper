#!/bin/zsh
export sim=$1 #change 
export sim_folder_name=$2 # change
export n_reps=$3 # change to 100 for full results
# export sequence=($4) # change to 1 for full results

for I in $(seq 1 ${n_reps})
do
  cd SimStudy/Results/${sim}/${sim_folder_name}/data 
  if [ ! -d "$I" ]; then
    mkdir $I
    cd $I
    for i in ${@[4,-1]} # change here 
    do
      mkdir res_nmtf_$i
      mkdir gfa_$i
      mkdir issvd_$i
      mkdir nmtf_$i
    done
  fi

  # move back into SimStudy
  cd ../../../../

  # generate data
  Rscript --vanilla data_gen.r  Results/${sim}/${sim_folder_name} $I

  # now analyse in R
  Rscript --vanilla methods_r.r  Results/${sim}/${sim_folder_name} $I

  # analyse in python
  python3 OtherMethods/methods_p.py Results/${sim}/${sim_folder_name} $I ${sim}

  # evaluate results in R
  Rscript --vanilla eval.r  Results/${sim}/${sim_folder_name} $I

  cd ..

done

# combine results
cd SimStudy/Results/${sim}/${sim_folder_name}
ls data/ > repeat_folders.txt
# move back into SimStudy
cd ../../../ 

#compile results
Rscript -e "rmarkdown::render('results.rmd')" Results/${sim}/${sim_folder_name}  Results/${sim}/${sim_folder_name}/repeat_folders.txt

mv results.pdf Results/${sim}/${sim_folder_name}