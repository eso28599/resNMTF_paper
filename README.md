# ResNMTF
 
This repository contains the files needed to recreate the results from the paper "Multi-view biclustering via non-negative matrix tri-factorisation".

It utilises our R packages `resnmtf` and `bisilhouette` which can be installed via:
```{r}
devtools::install_github("eso28599/bisilhouette") # requires devtools package to be installed
devtools::install_github("eso28599/resnmtf") 
```
If you want to use the ResNMTF rather than recreate the results, the `resnmtf` package is more appropriate, the repository of which can be found [here](https://github.com/eso28599/resnmtf).

The repository is split into three folders; 'SimStudy' which reproduces the simulation study results, 'RealData' which contains the and 'Exploration' which contains code to reproduce additional figures from the paper.

## SimStudy
This folder contains scripts to reproduce the simulation study results. It is set up for use on a HPC and has the following set up:

- Functions:
  - data_generation.r, synthetic data generation
  - evaluation_funcs.r, mutli-view biclustering evaluation
  - extra_funcs.r, useful functions
  - summary_funcs.r, functions to process results across repetitions
- Results, contains subdirectories set up to store data and results from simulations
- Scripts, contains scripts to run the simulation studies and extract results
- NoStability:
- OtherMethods:


For each study we have... 

For example, when we investigate increasing the number of biclusters for two views of data, the following folders and files are relevant.
- Results
    - bicl
       - bicl_2v
          - data, 
          - logs
          - sim_parameters.r, contains parameter settings for the simulation study
- Scripts: 
    - bicl_2v_script_one.sh, batch file repeating experiment 100 times
    - bicl_2v_script_two.sh, combines results and produces plots

In order to  ...
And the results (plots and .csv files) can then be found in Results/bicl/bicl_2v

### File names

- bicl, increases the number of biclusters present for data with
  - bicl_2v: 2 views
  - bicl_3v: 3 views
  - bicl_4v: 4 views
- exhaustive, 
  - no_overlap_3v5b,
  - overlap_3v5b
- indiv, increases the number of individuals present for data with
  - indiv_3v4b: 3 views and 4 biclusters
  - indiv_3v5b: 3 views and 5 biclusters
  - indiv_4v5b: 4 views and 5 biclusters

## RealData
 - contains datasets used in the paper, their processing and analysis.

## Exploration
- files to replicate convergence, correlation and JSD related figures.


# Citation
If you use our model in your work, please cite us with:

> Orme, E.S.C., Rodosthenous, T. and Evangelou, M., 2025. Multi-view biclustering via non-negative matrix tri-factorisation. arXiv preprint arXiv:2502.13698.

Or with the following bibtex entry:
```
@misc{resnmtf,
      title={Multi-view biclustering via non-negative matrix tri-factorisation}, 
      author={Ella S. C. Orme and Theodoulos Rodosthenous and Marina Evangelou},
      year={2025},
      eprint={2502.13698},
      archivePrefix={arXiv},
      primaryClass={stat.ME},
      url={https://arxiv.org/abs/2502.13698}, 
}
```