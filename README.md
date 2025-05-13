# ResNMTF
 
> [!info]
> To reproduce results from the paper "Multi-view biclustering via non-negative matrix tri-factorisation" (i.e. for the Bioinference conference reproductibility pirze) please run the `reproduce_results.rmd` file. 

This repository contains the files needed to recreate the results from the paper "Multi-view biclustering via non-negative matrix tri-factorisation".

It utilises our R packages `resnmtf` and `bisilhouette` which can be installed via:
```{r}
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
devtools::install_github("eso28599/bisilhouette")
devtools::install_github("eso28599/resnmtf") 
```
If you want to use the ResNMTF method rather than recreate the results, the `resnmtf` package may be more appropriate, the repository of which can be found [here](https://github.com/eso28599/resnmtf).

The repository is split into three folders; 'SimStudy' which reproduces the simulation study results, 'RealData' which contains the files to apply ResNMTF and other methods to several real datasets, and 'Exploration' which contains code to reproduce additional figures from the paper.

## SimStudy
This folder contains scripts to reproduce the simulation study results. It is set up for use on a HPC and has the following structure:
```
├── Functions:
│  ├── data_generation.r, functions to generate synthetic data 
│  ├── evaluation_funcs.r, functions to evaluate mutli-view biclustering
│  ├── extra_funcs.r, useful functions
│  ├── summary_funcs.r, functions to process results across repetitions
├── Results, contains subdirectories set up to store data and results from simulations
├── Scripts, contains scripts to run the simulation studies and extract results
├── NoStability, contains files to investigate the effect of stability analysis
├── OtherMethods, contains files to implement iSSVD, GFA and ResNMTF with no restrictions
├── data_gen.r, file to generate data
├── eval_ssvd_param.r, file to evaluate results from iSSVD parameter investigtion
├── eval.r, file to evaluate results
├── issvd_data_gen.r, file to generate data as described in iSSVD paper
├── issvd_results.rmd, file to produce results from iSSVD investigations
├── methods_r.r, file to apply methods in r
├── phi_results.rmd, file to produce results from phi investigations
├── results.rmd, file to produce results
```

Each study has corresponding scripts and folders to run and store the relevant data. For example, the following folders and files are specific to the "increasing biclusters over 2 views" study.
```
├── Results
│  ├── bicl
│  │  ├── bicl_2v
│  │  │  ├── data
│  │  │  │  ├── 001, synthetic data from the 1st repetition, true and found bicluster assignments
│  │  │  │  ├── ...
│  │  │  │  ├── 100, synthetic data from the 100th repetition, true and found bicluster assignments
│  │  │  ├── logs, log files for each repetition
│  │  │  ├── sim_parameters.r, contains parameter settings for the simulation study
├── Scripts: 
│  ├── bicl_2v_script_one.sh, batch file repeating experiment 100 times
│  ├── bicl_2v_script_two.sh, combines results and produces plots
```

Assuming the current directory is the ResNMTF_paper folder, the following lines of code are run on an HPC to run the "increasing biclusters over 2 views" simulation study: 
```
cd SimStudy/Scripts
qsub bicl_2v_script_one.sh
```
And once this has executed:
```
qsub bicl_2v_script_two.sh
```
And the results (plots and .csv files) can then be found in SimStudy/Results/bicl/bicl_2v. These include plots such as

### File names
The following subdirectories can be found in the SimStudy/Results folder: 
```
├── bicl, increases the number of biclusters present for data with
│  ├── bicl_2v: 2 views
│  ├── bicl_3v: 3 views
│  ├── bicl_4v: 4 views
├── exhaustive, increasing the non-exhaustivity rate 
│  ├── no_overlap_3v5b, 3 views and 5 biclusters with no overlap
│  ├── overlap_3v5b, 3 views and 5 biclusters with fixed overlap
├── indiv, increases the number of individuals present for data with
│  ├── indiv_3v4b: 3 views and 4 biclusters
│  ├── indiv_3v5b: 3 views and 5 biclusters
│  ├── indiv_4v5b: 4 views and 5 biclusters
├── issvd_data_gen, increases the number of individuals present for data 
│  ├── case1: ??
│  ├── case2: ??
├── issvd_param_mean100, investigates the effect of a hyperparameter on iSSVD performance for data with strong signal ($\mu=100$) and
│  ├── issvd_param3: 3 views and 3 biclusters
│  ├── issvd_param4: 3 views and 4 biclusters
│  ├── issvd_param: 3 views and 5 biclusters
├── issvd_param_signal5, investigates the effect of a hyperparameter on iSSVD performance for data with strong signal ($\mu=100$) and
│  ├── issvd_param3: 3 views and 3 biclusters
│  ├── issvd_param4: 3 views and 4 biclusters
│  ├── issvd_param5: 3 views and 5 biclusters
│  ├── issvd_param6: 3 views and 6 biclusters
├── noise, increases the level of noise present for data with
│  ├── noise_3v4b: 3 views and 4 biclusters
│  ├── noise_3v5b: 3 views and 5 biclusters
│  ├── noise_4v5b: 4 views and 5 biclusters
├── overlap, increasing the non-exhaustivity rate 
│  ├── exhaustive_3v5b, 3 views and 5 biclusters with no overlap
│  ├── non_exhaust_3v5b, 3 views and 5 biclusters with fixed overlap
├── signal, increases the level of signal present for data with
│  ├── signal3: 3 views and 4 biclusters
│  ├── signal4: 3 views and 5 biclusters
│  ├── signal4_all: 4 views and 5 biclusters
│  ├── signal5: 3 views and 5 biclusters
│  ├── signal6: 3 views and 5 biclusters
├── views, increases the number of views present for data with
│  ├── views_3b: 2 views
│  ├── views_4b: 3 views
│  ├── views_5b: 4 views
│  ├── views_5b_nostab: 4 views
```

## RealData
Contains datasets used in the paper, their processing and analysis.

## Exploration
Contains files to replicate convergence, correlation and JSD related figures.


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
