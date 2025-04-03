# ResNMTF
 
This repository contains the files needed to recreate the results from the paper ... . 

It utilises our R packages `resnmtf` and `bisilhouette` which can be installed via:
```{r}

```

The repository is split into three folders; 'SimStudy' which reproduces the simulation study results, 'RealData' which contains the and 'Exploration' which contains code to reproduce additional figures from the paper.

## SimStudy
This folder contains scripts to reproduce the simulation study results. It is set up for use on a HPC. 


- Functions:
  - data_generation.r, synthetic data generation
  - evaluation_funcs.r, mutli-view biclustering evaluation
  - extra_funcs.r, useful functions
  - summary_funcs.r, functions to process results across repetitions
- Results, contains subdirectories set up to store data and results from simulations
- Scripts, contains scripts to run the simulation studies and extract results
- NoStability:
- OtherMethods:


For each study. For example, when we investigate increasing the number of biclusters for two views of data
- Results
    - bicl
       - bicl_2v
          - data
          - logs
          - sim_parameters.r
- Scripts: 
    - bicl_2v_script_one.sh
    - bicl_2v_script_two.sh

And the results (plots and .csv files) can then be found in Results/bicl/bicl_2v




## Exploration
- files to replicate convergence, correlation and JSD results.

## RealData
 - contains datasets used in the paper, their processing and analysis.

