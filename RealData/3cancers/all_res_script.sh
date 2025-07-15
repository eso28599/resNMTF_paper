# python methods
conda activate resnmtf_env
python RealData/3cancers/3cancers_python_analysis.py
conda deactivate

# R based other methods
Rscript --vanilla RealData/other_results.r 3cancers
# combine results

# Rscript --vanilla RealData/combine_dis_metrics.r 3sources 3