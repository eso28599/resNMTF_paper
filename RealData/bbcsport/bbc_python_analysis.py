#3sources python analysis
from pandas import ExcelWriter
import pandas as pd
import numpy as np
from iSSVD.functions import issvd
from mvlearn.cluster import MultiviewCoRegSpectralClustering
from numpy.random import seed
import sys
sys.path.insert(0, '/rds/general/user/eso18/home/resNMTF_paper/SimStudy/OtherMethods')
from python_utils import save_csvs, fix_row_clusts, fix_col_clusts

file_path = "RealData/bbcsport"
data_views = pd.ExcelFile(file_path +  "/bbc_data_processed.xlsx")
data = [np.array(pd.read_excel(data_views, sheet)).transpose() for sheet in data_views.sheet_names]
n_views = len(data)
n_vars = [view.shape[1] for view in data]
n_samps = data[0].shape[0]

for i in np.arange(5):
        seed(10 + i)
        row_issvd_filename = "RealData/bbcsport/issvd_res/" +  str(i) + "_row_clusts"
        col_issvd_filename =  "RealData/bbcsport/issvd_res/" +  str(i) + "_col_clusts"
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
        save_csvs(row_clusts, row_issvd_filename)
        save_csvs(col_clusts, col_issvd_filename)
        
# now look at mvlearn 
for i in np.arange(5):
  seed = 10 + i
  row_mvclust_filename = "RealData/bbcsport/mvlearn_res/" + str(i) + "row_mvclusts"
  mv_spectral = MultiviewCoRegSpectralClustering(n_clusters=5,
      random_state=seed, n_init=100)
  mv_clusters = mv_spectral.fit_predict(data)
  num_clusters = mv_clusters.max() + 1
  row_cluster_mvc = np.zeros((mv_clusters.size, num_clusters), dtype=int)
  row_cluster_mvc[np.arange(mv_clusters.size), mv_clusters] = 1
  full_rows = [pd.DataFrame(row_cluster_mvc) for _ in np.arange(n_views)]
  save_csvs(full_rows, row_mvclust_filename)
