#!/usr/bin/env python
# coding: utf-8

# In[26]:


import scanpy as sc
import ucdeconvolve as ucd
import sys
import logging
import pandas as pd
from anndata import AnnData


# In[37]:


ucd.api.authenticate("uc_tA1ilU6JEWuUsOMUO8ARqowQyGY8XwJFqkJvn9chbtuNAlKS")


# In[27]:


path = "~/Documents/Spatial/Bustelo_OV/UCD_input/"
sample = sys.argv[1]


# In[29]:


#load counts exported from R and create annData
counts = pd.read_csv(path + sample + "_data_input.csv", index_col=0)
image = pd.read_csv(path + sample + "_coords_input.csv", index_col=0)


# In[30]:


counts.shape


# In[31]:


image.shape


# In[32]:


adata_vis = AnnData(counts.T)
adata_vis.var["SYMBOL"]=adata_vis.var_names
adata_vis.obs["array_row"]=image["row"]
adata_vis.obs["array_col"]=image["col"]
adata_vis.obsm["spatial"]=image[["imagerow", "imagecol"]].to_numpy()
adata_vis.obs["sample"]=image["tissue"]


# In[33]:


image[["imagerow", "imagecol"]].to_numpy().shape


# In[ ]:


ucd.tl.base(adata_vis,  verbosity = logging.DEBUG, split=True, propagate=True, use_raw=False)


# In[ ]:


out_path= "~/Documents/Spatial/Bustelo_OV/UCD_output/"


# In[35]:


predictions = ucd.utils.read_results(adata_vis, category="primary")
predictions.to_csv(out_path + sample + "_primary_propagated.csv", sep='\t')
predictions2 = ucd.utils.read_results(adata_vis, category="cancer")
predictions2.to_csv(out_path + sample + "_cancer.csv", sep='\t')

