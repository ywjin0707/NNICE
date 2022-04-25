"""
Created on Thu Apr 21 14:38:31 2022

@author: ywjin0707

scRNA-seq datasets
sample bulk dataset (TCGA) 'D:/OneDrive/210831_Bioinformatics/DATA/TCGA/TCGA_GDC_HTSeq_Counts.txt'
GRCh37.87 download link: http://ftp.ensembl.org/pub/grch37/release-100/gtf/homo_sapiens/
gtfpath = 'D:/OneDrive/210831_Bioinformatics/DATA/Gene Lists/Homo_sapiens.GRCh37.87.gtf'
directory = 'D:/OneDrive/210831_Bioinformatics/DATA/b_cells_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/'
filepath = 'D:/OneDrive/210831_Bioinformatics/DATA/TCGA/TCGA_GDC_HTSeq_Counts.txt'
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import pyranges as pr
import matplotlib.pyplot as plt
import seaborn as sns

#%%
# gtfpath = 'D:/OneDrive/210831_Bioinformatics/DATA/Gene Lists/Homo_sapiens.GRCh37.87.gtf'
gtfpath = 'C:/Users/yw_ji/OneDrive/210831_Bioinformatics/DATA/Gene Lists/Homo_sapiens.GRCh37.87.gtf'
# directory = 'D:/OneDrive/210831_Bioinformatics/DATA/b_cells_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/'
directory = 'C:/Users/yw_ji/210831_Bioinformatics/DATA/b_cells_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/'
# filepath = 'D:/OneDrive/210831_Bioinformatics/DATA/TCGA/TCGA_GDC_HTSeq_Counts.txt'
filepath = 'C:/Users/yw_ji/210831_Bioinformatics/DATA/TCGA/TCGA_GDC_HTSeq_Counts.txt'


bkdata = load_data(filepath)
bkdata = bkdata.T ### Need to fix in load_data
bkdata = generate_QC_metrics(bkdata, plot=False)
sns.jointplot(x='total_counts', y='n_genes_by_counts', height=8, data=bkdata.obs, kind='scatter')
sns.jointplot(x='total_counts', y='pct_counts_mito', height=8, data=bkdata.obs, kind='scatter')
sns.jointplot(x='total_counts', y='pct_counts_ribo', height=8, data=bkdata.obs, kind='scatter')
sns.jointplot(x='total_counts', y='pct_counts_mribo', height=8, data=bkdata.obs, kind='scatter')
bkdata = filter_by_QC_metrics(bkdata,  stds = 2, th_mito = 5, th_mribo = None, plot=False)
sns.jointplot(x='total_counts', y='n_genes_by_counts', height=8, data=bkdata.obs, kind='scatter')
sns.jointplot(x='total_counts', y='pct_counts_mito', height=8, data=bkdata.obs, kind='scatter')
sns.jointplot(x='total_counts', y='pct_counts_ribo', height=8, data=bkdata.obs, kind='scatter')
sns.jointplot(x='total_counts', y='pct_counts_mribo', height=8, data=bkdata.obs, kind='scatter')

scdata = load_data(directory)
scdata = generate_QC_metrics(scdata, plot=False)
sns.jointplot(x='total_counts', y='n_genes_by_counts', height=8, data=scdata.obs, kind='scatter')
sns.jointplot(x='total_counts', y='pct_counts_mito', height=8, data=scdata.obs, kind='scatter')
sns.jointplot(x='total_counts', y='pct_counts_ribo', height=8, data=scdata.obs, kind='scatter')
sns.jointplot(x='total_counts', y='pct_counts_mribo', height=8, data=scdata.obs, kind='scatter')
scdata = filter_by_QC_metrics(scdata,  stds = 2, th_mito = 5, th_mribo = None, plot=False)
sns.jointplot(x='total_counts', y='n_genes_by_counts', height=8, data=scdata.obs, kind='scatter')
sns.jointplot(x='total_counts', y='pct_counts_mito', height=8, data=scdata.obs, kind='scatter')
sns.jointplot(x='total_counts', y='pct_counts_ribo', height=8, data=scdata.obs, kind='scatter')
sns.jointplot(x='total_counts', y='pct_counts_mribo', height=8, data=scdata.obs, kind='scatter')

gtf = pd.read_table(gtfpath, header=None, comment='#', dtype=object)
col8 = gtf[8].str.split(r'\s')
#%%

def load_gtf(path):
    dtypes = {
        "Chromosome": "category",
        "Feature": "category",
        "Strand": "category"
    }
    
    names = "Chromosome Source Feature Start End Score Strand Frame Attribute".split()
    
    gtf = pd.read_csv(
        path,
        sep='\t',
        header=None,
        comment='#',
        names=names,
        dtype=dtypes)
    return gtf
    
def process_gtf(gtf):
    # gtf = gtf.loc[gtf['Chromosome'].isin([str(i) for i in range(1,23)] + ['X','Y','MT'])]

    rowdicts = []
    for l in gtf.Attribute:
        rowdict = {}
        # l[:-1] removes final ";" cheaply
        for k, v in (kv.replace('"', '').split(None, 1) for kv in l[:-1].split("; ")):
            if k in ['gene_id', 'gene_name']:
                if k not in rowdict:
                    rowdict[k] = v
                elif k in rowdict and isinstance(rowdict[k], list):
                    rowdict[k].append(v)
                else:
                    rowdict[k] = [rowdict[k], v]

        rowdicts.append({
            k: ','.join(v) if isinstance(v, list) else v
            for k, v in rowdict.items()
        })
    extra = pd.DataFrame.from_dict(rowdicts).set_index(gtf.index)
    gtf = gtf.drop(["Attribute"], axis=1)
    df = pd.concat([gtf, extra], axis=1, sort=False)
    # df.loc[:, "Start"] = df.Start - 1
    return df

        
def load_data(path, sampleID=None, celltype=None):
    """
    Loads data types including:
        - 10x Genomics scRNA-seq data from directory. 
        - .csv file format
        - .txt file format
    Annotate sample and celltype information. 
    
    Parameters
    ----------
    directory : str
        Relative or absolute path of directory with gene expression data. (ex. '../DATA/b_cell/filtered_matrices_mex/hg19/')
    sampleID : str, int, optional
        Sample identification, if known. Default is None. 
    celltype : str, optional
        Cell type identity, if known. Default is None. 

    Returns
    -------
    mydata : AnnData class object
        Count matrix with sampleID and celltype annotations.

    """
    if os.path.isdir(path):
        mydata = sc.read_10x_mtx(path)
    else:
        try:
            mydata = ad.read_csv(path)
        except ValueError:
            mydata = ad.read_text(path)
    mydata.uns['filepath'] = path
    if sampleID:
        mydata.obs['sampleID'] = str(sampleID)
    if celltype:
        mydata.obs['celltype'] = celltype
        mydata.obs.set_index(pd.Index([celltype+'_'+rn for rn in mydata.obs.index]), inplace=True)
    return mydata




def generate_QC_plots(mydata, filtered=False):
    """
    Make 3 joint grid plots showing quality control metrics for cells/samples:
        - Total counts X Number of features
        - Total counts X Percent mitochondrial genes
        - Total counts X Percent mitochondrial ribosomal genes

    Parameters
    ----------
    data : AnnData class object
        Count matrix with annotations for QC metrics.
    filtered : bool, optional
        Whether the count matrix was filtered already. The default is False.

    Returns
    -------
    None.

    """
    filepath = mydata.uns['filepath']
    if filtered:
        suffix = '_filtered.png'
    else:
        suffix = '.png'
    sns.jointplot(x='total_counts', y='n_genes_by_counts', 
                  height=8, data=mydata.obs, kind='scatter').savefig(filepath + 'totalcounts_totalgenes' + suffix)
    sns.jointplot(x='total_counts', y='pct_counts_mito', 
                  height=8, data=mydata.obs, kind='scatter').savefig(filepath + 'totalcounts_pctmito' + suffix)
    sns.jointplot(x='total_counts', y='pct_counts_ribo', 
                  height=8, data=mydata.obs, kind='scatter').savefig(filepath + 'totalcounts_pctribo' + suffix)
    sns.jointplot(x='total_counts', y='pct_counts_mribo', 
                  height=8, data=mydata.obs, kind='scatter').savefig(filepath + 'totalcounts_pctmribo' + suffix)


def generate_QC_metrics(mydata, plot=True):
    """
    Annotate mitochondrial (mito), ribosomal (ribo), and mitochondrial-ribosomal (mribo) genes.
    Calculate QC metrics - refer to scanpy docs for more details (https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.calculate_qc_metrics.html#scanpy.pp.calculate_qc_metrics). 
    Then plots QC metrics (optional). 
                                                                  
    Parameters
    ----------
    mydata : AnnData class object
        Loaded gene expression data.
    plot : bool, optional
        Whether to plot QC metrics in the original input data directory. The default is True.

    Returns
    -------
    mydata : AnnData class object
        Count matrix with QC metrics.

    """
    mydata.var['mito'] = mydata.var.index.str.match('^MT-')
    mydata.var['ribo'] = mydata.var.index.str.startswith(('RPL','RPS'))
    mydata.var['mribo'] = mydata.var.index.str.startswith(('MRPL','MRPS'))
    sc.pp.calculate_qc_metrics(mydata, qc_vars=['mito','ribo','mribo'], inplace=True)
    if plot:
        generate_QC_plots(mydata, filtered=False)
    return mydata


def filter_by_QC_metrics(mydata,  stds = 2, th_mito = 5, th_mribo = None, plot=True):
    """
    Filters cells by 4 metrics: 
        - total read counts for each cell/sample within 2 (default) standard deviations from the mean across all cells/samples;
        - total number of non-zero features for each cell/sample within 2 (default) standard deviations from the mean across all cells/samples;
        - percentage of total transcripts from mitochondrial genes; and
        - percentage of total transcripts from mitochondrial ribosomal genes. 

    Parameters
    ----------
    mydata : AnnData class object
        AnnData with annotations for QC metrics.
    stds : float
        Number of standard deviations around the mean metric value to keep. The default is 2. 
    th_mito : float
        Threshold for mitochondrial gene content from 0-100. The default is 5.
    th_mribo : float, optional
        Threshold for mitochondrial ribosomal gene content from 0-100. The default is None.

    Returns
    -------
    mydata : AnnData class object
        Filtered scRNA-seq count matrix.

    """
    # sc.pp.filter_cells(mydata, **filter_cells_args) # inplace=True by default
    # sc.pp.filter_genes(mydata, **filter_genes_args) # inplace=True by default
    mean_total_counts = mydata.obs['total_counts'].mean()
    std_total_counts = mydata.obs['total_counts'].std()
    mean_n_genes_by_counts = mydata.obs['n_genes_by_counts'].mean()
    std_n_genes_by_counts = mydata.obs['n_genes_by_counts'].std()
    mydata = mydata[(mydata.obs.total_counts > (mean_total_counts - std_total_counts*stds)) & (mydata.obs.total_counts < (mean_total_counts + std_total_counts*stds)), :]
    mydata = mydata[(mydata.obs.n_genes_by_counts > (mean_n_genes_by_counts - std_n_genes_by_counts*stds)) & (mydata.obs.n_genes_by_counts < (mean_n_genes_by_counts + std_n_genes_by_counts*stds)), :]
    mydata = mydata[mydata.obs.pct_counts_mito < th_mito, :]
    if th_mribo:
        mydata = mydata[mydata.obs.pct_counts_mribo < th_mribo, :]
    if plot:
        generate_QC_plots(mydata, qc=True)
    return mydata
    

def concat(*args):
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
    res = ad.concat([mydata for mydata in args])
    return res

