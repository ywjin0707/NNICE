import tensorflow as tf
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from model import MinMaxNorm

def FeatureList(paths: list) -> list:
    features = None
    for path in paths:
        mydata = pd.read_csv(path, index_col = 0)
        if features == None:
            features = set(mydata.index.values.tolist())
        else:
            features = features.intersection(set(mydata.index.values.tolist()))
    features = list(features)
    features.sort()
    return features


class DataPreprocess():

    def __init__(self, datadir, celltypes, bkdata_path, features):
        '''
        Creates preprocessed instance of input data
        scdata should be in matrix.mtx within specified folders along with barcodes.tsv and genes.tsv
        bkdata should have sample names as columns and gene names as rows
        gene_list should have no row or column names/index
        '''
        self.datadir = datadir
        self.celltypes = celltypes
        self.scdata = self.load_scdata(self.datadir, self.celltypes)
        self.bkdata = pd.read_csv(bkdata_path, index_col=0)
        # If there is input gene list, filter out genes not in bkdata or scdata
        if features is None:
            self.features = self.bkdata.index.drop_duplicates()
        else:
            self.features = features
        # Filter out genes not in gene list
        self.scdata = self.scdata[:,self.scdata.var_names.isin(self.features)]
        sc.pp.normalize_total(self.scdata, target_sum=1e6) # normalize to sum to 1,000,000
        # sc.pp.regress_out(scdata, ['total_counts'], n_jobs=1)
        # Transpose, filter out genes not in gene list, then sort column (by gene name)
        self.bkdata = self.bkdata.T
        self.bkdata = self.bkdata.loc[:,self.bkdata.columns.isin(self.features)].sort_index(axis=1)
        self.bkdata = self.bkdata.values.astype(float)

    def load_scdata(self, data_directories, cell_types):
        # Read and merge 10X Genomics scRNA-seq data
        scdata = None
        print('Loading single cell dataset')
        for d, c in zip(tqdm(data_directories), cell_types):
            x = sc.read_10x_mtx(d)
            x.obs['celltype'] = [c]*len(x.obs.index)
            # Change each observation (cell) name to celltype + barcode
            x.obs.set_index(pd.Index([c+'_'+rn[:-2] for rn in x.obs.index]), inplace=True)
            if scdata is not None:
                scdata = ad.concat([scdata, x])
            else:
                scdata = x
        # Filter out cells and genes
        sc.pp.filter_cells(scdata, min_genes=200)
        sc.pp.filter_genes(scdata, min_cells=1)
        # Search for prefix "MT-" (mitochondrial genes) and make new column in variable annotations
        # Search for prefix "RPL/RPS" for ribosomal genes and "MRPL/MRPS" for mitochondrial ribosomal genes
        scdata.var['mito'] = scdata.var.index.str.match('^MT-')
        scdata.var['ribo'] = scdata.var.index.str.startswith(('RPL','RPS'))
        scdata.var['mribo'] = scdata.var.index.str.startswith(('MRPL','MRPS'))
        # Calculate QC metrics as per McCarthy et al., 2017 (Scater)
        sc.pp.calculate_qc_metrics(scdata, qc_vars=['mito','ribo', 'mribo'], inplace=True)
        # Plot QC metrics
        # sns.jointplot(x='total_counts', y='n_genes_by_counts', height=8, data=scdata.obs,
        #     kind='scatter', hue='celltype')
        # sns.jointplot(x='total_counts', y='pct_counts_mito', height=8, data=scdata.obs,
        #     kind='scatter', hue='celltype')
        # sns.jointplot(x='total_counts', y='pct_counts_ribo', height=8, data=scdata.obs,
        #     kind='scatter', hue='celltype')
        # sns.jointplot(x='total_counts', y='pct_counts_mribo', height=8, data=scdata.obs,
        #     kind='scatter', hue='celltype')
        # plt.show()
        # Filter out cells with >5% of counts from mitochondria and mitoribosome
        # scdata = scdata[scdata.obs.pct_counts_ribo > 30, :]
        scdata = scdata[scdata.obs.pct_counts_mito < 5, :]
        scdata = scdata[scdata.obs.pct_counts_mribo < 1, :]
        return scdata

    def __call__(self, whichdata, batch_size=1):
        if whichdata == 'scdata':
            out = []
            print('Dividing single cell dataset into cell types')
            for c in tqdm(self.celltypes):
                scdata_ = self.scdata[self.scdata.obs.celltype==c].to_df().sort_index(axis=1)
                # Add to row index 0 a cell with no gene expression (all zeros)
                # zeros = pd.DataFrame(np.zeros((1,scdata_.shape[1])), columns=scdata_.columns.values)
                # Expand into batch dimension and repeat 2-D tensor by # of samples per mini batch
                # scdata_ = tf.tile(tf.expand_dims(pd.concat([zeros,scdata_]), axis=0), [batch_size,1,1])
                out.append(scdata_)
        elif whichdata == 'bkdata':
            out = tf.data.Dataset.from_tensor_slices(MinMaxNorm(tf.math.log1p(self.bkdata)))
            out = out.shuffle(100).batch(1)
        elif whichdata == 'genelist':
            out = self.features
        else:
            raise ValueError('Choose only one of the following: "scdata", "bkdata", or "genelist"')
        return out
