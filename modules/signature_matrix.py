# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 15:39:55 2021

@author: yw_ji
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
# import seaborn as sns
# import matplotlib.pyplot as plt
from tqdm import tqdm


os.chdir('C:/Users/yw_ji/Documents/MSc Thesis/dice/')

#%%
wu_meta = pd.read_csv('../DATA/TNBC/Wu_EMBO_metadata.csv', skiprows= lambda t: t in [1])
# wu_meta.celltype_final.unique()

wu_meta.loc[wu_meta.celltype_final == 'iCAFs', 'celltype_final'] = 'CAF'
wu_meta.loc[wu_meta.celltype_final == 'myCAFs', 'celltype_final'] = 'CAF'
wu_meta.loc[wu_meta.celltype_final == 'Plasma_Cells', 'celltype_final'] = 'B_cell'
wu_meta.loc[wu_meta.celltype_final == 'Epithelial_Basal', 'celltype_final'] = 'Cancer'
wu_meta.loc[wu_meta.celltype_final == 'T_cells_unassigned', 'celltype_final'] = 'OMIT'
wu_meta.loc[wu_meta.celltype_final == 'CD8+ T-cells', 'celltype_final'] = 'CD8+_T_cell'
wu_meta.loc[wu_meta.celltype_final == 'NKT cells', 'celltype_final'] = 'NK_cell'
wu_meta.loc[wu_meta.celltype_final == 'NK cells', 'celltype_final'] = 'NK_cell'
wu_meta.loc[wu_meta.celltype_final == 'T-cells Cycling', 'celltype_final'] = 'OMIT'
wu_meta.loc[wu_meta.celltype_final == 'CD4+ T-cells', 'celltype_final'] = 'CD4+_T_cell'
wu_meta.loc[wu_meta.celltype_final == 'T-Regs', 'celltype_final'] = 'Treg'
wu_meta.loc[wu_meta.celltype_final == 'Tfh cells', 'celltype_final'] = 'CD4+_T_cell'
wu_meta.loc[wu_meta.celltype_final == 'dPVL', 'celltype_final'] = 'PVL'
wu_meta.loc[wu_meta.celltype_final == 'imPVL', 'celltype_final'] = 'PVL'
wu_meta.loc[wu_meta.celltype_final == 'B_Cells', 'celltype_final'] = 'B_cell'
wu_meta.loc[wu_meta.celltype_final == 'Myoepithelial', 'celltype_final'] = 'Epithelial'
wu_meta.loc[wu_meta.celltype_final == 'Epithelial_Basal_Cycling', 'celltype_final'] = 'Cancer'
wu_meta.loc[wu_meta.celltype_final == 'Epithelial_Luminal_Mature', 'celltype_final'] = 'Epithelial'

wu_meta = wu_meta.set_index('NAME')
# del wu_meta.index.name

#%%
scdata = sc.read_10x_mtx('C:/Users/yw_ji/Documents/MSc Thesis/DATA/TNBC/counts_matrix/')
scdata.obs = pd.DataFrame(scdata.obs.merge(wu_meta, left_index=True, right_index=True).loc[:,['celltype_final','patientID']])
scdata.obs.columns = ['celltype', 'patientID']
scdata.obs.set_index(pd.Index([c+'_'+rn[-16:-1] for c, rn in zip(scdata.obs.celltype, scdata.obs.index)]), inplace=True)

#%%
data_directories = ["../DATA/b_cells_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "../DATA/cd4_t_helper_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "../DATA/cd14_monocytes_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "../DATA/cd34_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "../DATA/cd56_nk_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "../DATA/cytotoxic_t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "../DATA/memory_t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "../DATA/naive_cytotoxic_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "../DATA/naive_t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "../DATA/regulatory_t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/"]
cell_types = ['B_cell','CD4+_T_cell','Myeloid','Myeloid','NK_cell','CD8+_T_cell','CD4+_T_cell','CD8+_T_cell','CD4+_T_cell','Treg']

#%%
data_paths = ['../DATA/TCGA/TCGA_GDC_HTSeq_TPM.csv',
              '../DATA/METABRIC/METABRIC.csv',
              '../DATA/SDY67/SDY67_468.csv']

#%%
def FeatureList(paths: list) -> list:
    features = None
    for path in tqdm(paths):
        mydata = pd.read_csv(path, index_col = 0)
        if features == None:
            features = set(mydata.index.values.tolist())
        else:
            features = features.intersection(set(mydata.index.values.tolist()))
    features = list(features)
    features.sort()
    return features

#%%
features = FeatureList(data_paths)

#%%
for d, c in zip(tqdm(data_directories), cell_types):
    x = sc.read_10x_mtx(d)
    x.obs['celltype'] = [c]*len(x.obs.index)
    x.obs['patientID'] = ['P0']*len(x.obs.index)
    # Change each observation (cell) name to celltype + barcode
    x.obs.set_index(pd.Index([c+'_'+rn[:-2] for rn in x.obs.index]), inplace=True)
    scdata = ad.concat([scdata, x])
    
# dim(scdata) = (118926, 25584)

#%%
# Filter out cells and genes
sc.pp.filter_cells(scdata, min_genes=200)
sc.pp.filter_genes(scdata, min_cells=1)
# Search for prefix "MT-" (mitochondrial genes) and make new column in variable annotations
# Search for prefix "RPL/RPS" for ribosomal genes and "MRPL/MRPS" for mitochondrial ribosomal genes
scdata.var['mito'] = scdata.var.index.str.match('^MT-')
scdata.var['ribo'] = scdata.var.index.str.startswith(('RPL','RPS'))
scdata.var['mribo'] = scdata.var.index.str.startswith(('MRPL','MRPS'))
# Calculate QC metrics as per McCarthy et al. 2017 (Scater)
sc.pp.calculate_qc_metrics(scdata, qc_vars=['mito','ribo', 'mribo'], inplace=True)

# dim(scdata) = (118859, 24690)

#%%
# Filter out cells with >5% of counts from mitochondria and mitoribosome
# scdata = scdata[scdata.obs.pct_counts_ribo > 30, :]
scdata = scdata[scdata.obs.pct_counts_mito < 5, :]
scdata = scdata[scdata.obs.pct_counts_mribo < 1, :]

# dim(scdata) = (105733, 24690)

#%%
# Filter cells with abstract cell type label
scdata = scdata[scdata.obs['celltype'] != 'OMIT']
# Filter genes not in any of the bulk GEPs
scdata = scdata[:,scdata.var_names.isin(features)]

# dim(scdata) = (104905, 14444)

#%%
# Curate mean expression levels of genes over all single cell samples in each cell type separately
mean_counts = None
for c in scdata.obs.celltype.unique():
    tmp = scdata[scdata.obs['celltype']==c]
    sc.pp.calculate_qc_metrics(tmp, inplace=True)
    if mean_counts is not None:
        mean_counts = pd.concat([mean_counts, tmp.var['mean_counts']], axis=1)
    else:
        mean_counts = tmp.var['mean_counts']

#%%
bins = pd.DataFrame({'max':np.max(mean_counts, axis=1)})
bins['binsize'] = bins['max']/3

#%%
tmp = scdata.X.toarray()
for i in range(len(scdata.obs)):
    tmp[i,:] = np.where(tmp[i,:] < bins['binsize'], 0, 
                        np.where(tmp[i,:] >= bins['max'], 2, 1))

#%%
scdata.X = tmp
del tmp

#%%
sc.pl.heatmap(scdata, var_name=scdata.var)