#iSSVD application  
import sys
path_to_sim_folder = str(sys.argv[1])
batch_folder = str(sys.argv[2])
from iSSVD.functions import issvd
import numpy as np
import pandas as pd
# import OtherMethods.python_utils
from OtherMethods.python_utils import fix_row_clusts, fix_col_clusts, save_xls

#method_idx = list(np.round(np.linspace(0.1, 1, 10), 2)) #only thing to change
method_idx = list(np.arange(1,11))

data_name = path_to_sim_folder + "/data/" + batch_folder 

for j in method_idx:
    data_views = pd.ExcelFile(data_name + "/res_nmtf_" + str(j) + "/data.xlsx")
    data = [np.array(pd.read_excel(data_views, sheet)) for sheet in data_views.sheet_names]
    #save row clustering
    n_views = len(data)
    n_vars = [view.shape[1] for view in data]
    n_samps = data[0].shape[0]
    row_issvd_filename = data_name + "/issvd_" + str(j) + "/row_clusts.xlsx"
    col_issvd_filename = data_name + "/issvd_" + str(j) + "/col_clusts.xlsx"
    iSSVD_applied = issvd(data, standr=False, pointwise=True,steps=100,size=0.5,
                vthr = 0.7,ssthr=[0.6,0.65],nbicluster=10,rows_nc=True,cols_nc=True,col_overlap=True
                ,row_overlap=True,pceru=0.1,pcerv=0.1,merr=0.0001,iters=100)           
    n_clusts = iSSVD_applied['N']
    if n_clusts == 0:
        row_clusts = [pd.DataFrame([0 for i in np.arange(k)]) for k in [n_samps for j in np.arange(n_views)]]
        col_clusts = [pd.DataFrame([0 for i in np.arange(k)]) for k in n_vars]
    else:
        row_clusts = fix_row_clusts(iSSVD_applied['Sample_index'], n_views, n_samps)
        col_clusts = fix_col_clusts(iSSVD_applied['Variable_index'], n_views, n_clusts, n_vars)
    save_xls(row_clusts, row_issvd_filename)
    save_xls(col_clusts, col_issvd_filename)
