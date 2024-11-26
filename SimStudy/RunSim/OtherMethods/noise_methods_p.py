#iSSVD application  
import sys
path_to_sim_folder = str(sys.argv[1])
batch_folder = str(sys.argv[2])
import pandas as pd
import numpy as np
from iSSVD.functions import issvd
from pandas import ExcelWriter

method_idx = list(np.arange(1,11)) + list(np.arange(20,110,10))#only thing to change

def save_xls(list_dfs, xls_path):
    with ExcelWriter(xls_path) as writer:
        for n, df in enumerate(list_dfs):
            df.to_excel(writer,'sheet%s' % n, header=True, index= False)
    
def cluster_assignment(clust_membership, n_var):
    #eturn [pd.DataFrame([np.arange(n_var) in clust for clust in clust_membership])]
    return [[int(i in clust) for clust in clust_membership]  for i in np.arange(n_var)]

def fix_col_clusts(clust_membership, n_views, n_clust, n_var):
    #clust_membership[i][j] is the i^th cluster in the j^th row - so the i^th column in the j^th view
    new_clust = [[clust_membership[i][j] for i in np.arange(n_clust)] for j in np.arange(n_views)]
    #in new_clust - the length of the list is the n_views, 
    #item in zip(new_clust, n_var) is the set of e.g. 10 clusters for the ith view
    #and the number of variables we are clustering here
    return [pd.DataFrame(cluster_assignment(clust,n)) for (clust, n) in zip(new_clust, n_var)]

def fix_row_clusts(clust_membership, n_views, n_samp):
    new_clust = [cluster_assignment(clust_membership, n_samp) for i in np.arange(n_views)]
    return [pd.DataFrame(results) for results in new_clust]

data_name = path_to_sim_folder + "/data/" + batch_folder 
#repeat 5 times,
for j in method_idx:
    data_views = pd.ExcelFile(data_name + "/res_nmtf_" + str(j) + "/data.xlsx")
    data = [np.array(pd.read_excel(data_views, sheet)) for sheet in data_views.sheet_names]

    #save row clustering
    n_views = len(data)
    n_vars = [view.shape[1] for view in data]
    n_samps = data[0].shape[0]
    for i in np.arange(5):
        row_issvd_filename = data_name + "/issvd_" + str(j) + "/" +  str(i) + "_row_clusts.xlsx"
        col_issvd_filename = data_name + "/issvd_" + str(j) + "/" +  str(i) + "_col_clusts.xlsx"
        iSSVD_applied = issvd(data, standr=False,pointwise=True,steps=100,size=0.6,
                    vthr = 0.7,ssthr=[0.6,0.8],nbicluster=10,rows_nc=True,cols_nc=True,col_overlap=False
                    ,row_overlap=False,pceru=0.15,pcerv=0.16,merr=0.0001,iters=100)           
        n_clusts = iSSVD_applied['N']
        if n_clusts == 0:
            row_clusts = [pd.DataFrame([0 for i in np.arange(k)]) for k in [n_samps for j in np.arange(n_views)]]
            col_clusts = [pd.DataFrame([0 for i in np.arange(k)]) for k in n_vars]
        else:
            row_clusts = fix_row_clusts(iSSVD_applied['Sample_index'], n_views, n_samps)
            col_clusts = fix_col_clusts(iSSVD_applied['Variable_index'], n_views, n_clusts, n_vars)
        save_xls(row_clusts, row_issvd_filename)
        save_xls(col_clusts, col_issvd_filename)

