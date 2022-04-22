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
    if sampleID:
        scdata.obs['sampleID'] = str(sampleID)
    if celltype:
        scdata.obs['celltype'] = celltype
        scdata.obs.set_index(pd.Index([celltype+'_'+rn for rn in scdata.obs.index]), inplace=True)
    return scdata


def generate_QC_plots(scdata, filtered=False):
    """
    Make 3 joint grid plots showing quality control metrics for cells:
        - Total counts X Number of features
        - Total counts X Percent mitochondrial genes
        - Total counts X Percent mitochondrial ribosomal genes

    Parameters
    ----------
    scdata : AnnData class object
        scRNA-seq count matrix with annotations for QC metrics.
    filtered : bool, optional
        Whether the count matrix was filtered already. The default is False.

    Returns
    -------
    None.

    """
    filepath = scdata.uns['directory']
    if filtered:
        suffix = '_filtered.png'
    else:
        suffix = '.png'
    sns.jointplot(x='total_counts', y='n_genes_by_counts', 
                  height=8, data=scdata.obs, kind='scatter').savefig(filepath + 'totalcounts_totalgenes' + suffix)
    sns.jointplot(x='total_counts', y='pct_counts_mito', 
                  height=8, data=scdata.obs, kind='scatter').savefig(filepath + 'totalcounts_pctmito' + suffix)
    sns.jointplot(x='total_counts', y='pct_counts_ribo', 
                  height=8, data=scdata.obs, kind='scatter').savefig(filepath + 'totalcounts_pctribo' + suffix)
    sns.jointplot(x='total_counts', y='pct_counts_mribo', 
                  height=8, data=scdata.obs, kind='scatter').savefig(filepath + 'totalcounts_pctmribo' + suffix)


def generate_QC_metrics(scdata, plot=True):
    """
    Annotate mitochondrial (mito), ribosomal (ribo), and mitochondrial-ribosomal (mribo) genes.
    Calculate QC metrics - refer to scanpy docs for more details (https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.calculate_qc_metrics.html#scanpy.pp.calculate_qc_metrics). 
    Then plots QC metrics (optional). 
                                                                  
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
    scdata.var['mito'] = scdata.var.index.str.match('^MT-')
    scdata.var['ribo'] = scdata.var.index.str.startswith(('RPL','RPS'))
    scdata.var['mribo'] = scdata.var.index.str.startswith(('MRPL','MRPS'))
    sc.pp.calculate_qc_metrics(scdata, qc_vars=['mito','ribo','mribo'], inplace=True)
    if plot:
        generate_QC_plots(scdata, filtered=False)
    return scdata


def filter_by_QC_metrics(scdata,  stds = 2, th_mito = 5, th_mribo = 1, plot=True):
    """
    Filters cells by 4 metrics: 
        - total read counts for each cell within 2 (default) standard deviations from the mean across all cells;
        - total number of non-zero features for each cell within 2 (default) standard deviations from the mean across all cells;
        - percentage of total transcripts from mitochondrial genes; and
        - percentage of total transcripts from mitochondrial ribosomal genes. 

    Parameters
    ----------
    scdata : AnnData class object
        AnnData with annotations for QC metrics.
    stds : float
        Number of standard deviations around the mean metric value to keep. The default is 2. 
    th_mito : float
        Threshold for mitochondrial gene content from 0-100. The default is 5.
    th_mribo : float
        Threshold for mitochondrial ribosomal gene content from 0-100. The default is 1.

    Returns
    -------
    scdata : AnnData class object
        Filtered scRNA-seq count matrix.

    """
    # sc.pp.filter_cells(scdata, **filter_cells_args) # inplace=True by default
    # sc.pp.filter_genes(scdata, **filter_genes_args) # inplace=True by default
    mean_total_counts = scdata.obs['total_counts'].mean()
    std_total_counts = scdata.obs['total_counts'].std()
    mean_n_genes_by_counts = scdata.obs['n_genes_by_counts'].mean()
    std_n_genes_by_counts = scdata.obs['n_genes_by_counts'].std()
    scdata = scdata[(scdata.obs.total_counts > (mean_total_counts - std_total_counts*stds)) & (scdata.obs.total_counts < (mean_total_counts + std_total_counts*stds)), :]
    scdata = scdata[(scdata.obs.n_genes_by_counts > (mean_n_genes_by_counts - std_n_genes_by_counts*stds)) & (scdata.obs.n_genes_by_counts < (mean_n_genes_by_counts + std_n_genes_by_counts*stds)), :]
    scdata = scdata[scdata.obs.pct_counts_mito < th_mito, :]
    scdata = scdata[scdata.obs.pct_counts_mribo < th_mribo, :]
    if plot:
        generate_QC_plots(scdata, qc=True)
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