"""
Created on Thu Apr 21 14:38:31 2022

@author: ywjin0707
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns


def load_scdata(directory, sampleID, celltype):
    """
    Loads 10x Genomics scRNA-seq data from directory. 
    Annotate mitochondrial (mito), ribosomal (ribo), and mitochondrial-ribosomal (mribo) genes.
    Annotate sample and celltype information. 
    
    Parameters
    ----------
    directory : str
        Relative path of directory with scRNA-seq data. (ex. '../DATA/b_cell/filtered_matrices_mex/hg19/')
    sampleID : str, int, optional
        Sample identification, if known.
    celltype : str, optional
        Cell type identity, if known.

    Returns
    -------
    scdata : AnnData class object
        scRNA-seq count matrix with sampleID and celltype annotations.

    """
    scdata = sc.read_10x_mtx(directory)
    scdata.uns['directory'] = directory
    scdata.var['mito'] = scdata.var.index.str.match('^MT-')
    scdata.var['ribo'] = scdata.var.index.str.startswith(('RPL','RPS'))
    scdata.var['mribo'] = scdata.var.index.str.startswith(('MRPL','MRPS'))
    if sampleID:
        scdata.obs['sampleID'] = str(sampleID)
    if celltype:
        scdata.obs['celltype'] = celltype
        scdata.obs.set_index(pd.Index([celltype+'_'+rn for rn in scdata.obs.index]), inplace=True)
    return scdata

def generate_QC_metrics(scdata, plot=True):
    """
    

    Parameters
    ----------
    scdata : AnnData class object
        Loaded scRNA-seq data.
    plot : bool, optional
        Whether to plot QC metrics in the original input data directory. The default is True.

    Returns
    -------
    scdata : AnnData class object
        scRNA-seq count matrix with QC metrics.

    """
    sc.pp.calculate_qc_metrics(scdata, qc_vars=['mito','ribo','mribo'], inplace=True)
    if plot:
        filepath = scdata.uns['directory']
        sns.jointplot(x='total_counts', y='n_genes_by_counts', 
                      height=8, data=scdata.obs, kind='scatter', hue='celltype').savefig(filepath + 'totalcounts_totalgenes.png')
        sns.jointplot(x='total_counts', y='pct_counts_mito', 
                      height=8, data=scdata.obs, kind='scatter', hue='celltype').savefig(filepath + 'totalcounts_pctmito.png')
        sns.jointplot(x='total_counts', y='pct_counts_ribo', 
                      height=8, data=scdata.obs, kind='scatter', hue='celltype').savefig(filepath + 'totalcounts_pctribo.png')
        sns.jointplot(x='total_counts', y='pct_counts_mribo', 
                      height=8, data=scdata.obs, kind='scatter', hue='celltype').savefig(filepath + 'totalcounts_pctmribo.png')
    return scdata

def filter_by_QC_metrics(scdata, filter_cells_args = {'min_genes':500}, num_std = 2):
    """
    

    Parameters
    ----------
    scdata : TYPE
        DESCRIPTION.
    filter_cells_args : dict, optional
        Keyworded arguments for sc.pp.filter_cells function. The default is {'min_genes':500}.
    filter_genes_args : dict, optional *REMOVED*
        Keyworded arguments for sc.pp.filter_genes function. The default is {'min_cells':3}.

    Returns
    -------
    scdata : AnnData class object
        DESCRIPTION.

    """
    sc.pp.filter_cells(scdata, **filter_cells_args) # inplace=True by default
    # sc.pp.filter_genes(scdata, **filter_genes_args) # inplace=True by default
    scdata = scdata[scdata.obs.pct_counts_mito < 5, :]
    scdata = scdata[scdata.obs.pct_counts_mribo < 1, :]
    return scdata
    

def concat(x, y, z):
    """
    

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.

    Returns
    -------
    res : TYPE
        DESCRIPTION.

    """
    res = pd.concat([x,y,z])
    return res