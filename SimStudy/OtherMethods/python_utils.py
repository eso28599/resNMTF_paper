import numpy as np
import pandas as pd
from pandas import ExcelWriter


# python utils
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