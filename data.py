import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns


class DataPreprocess():
    def __init__(self, datadir, celltypes, bkdata_path, gene_list_path):
        '''
        Creates preprocessed instance of input data
        scdata should be in matrix.mtx within specified folders along with barcodes.tsv and genes.tsv
        bkdata should have sample names as columns and gene names as rows
        gene_list should have no row or column names/index
        '''
        self.datadir = datadir
        self.celltypes = celltypes
        self.scdata = load_scdata(self.datadir, self.celltypes)
        self.bkdata = pd.read_csv(bkdata_path, sep='\t')
        # If there is input gene list, filter out genes not in bkdata or scdata
        if gene_list_path is not None:
            self.genelist = pd.read_csv(gene_list_path, squeeze=True, names=['name'])
            self.genelist = self.genelist.drop_duplicates()
            self.genelist = self.genelist[self.genelist.isin(self.bkdata.index.drop_duplicates())]
        else:
            self.genelist = self.bkdata.index.drop_duplicates()
        self.genelist = self.genelist[self.genelist.isin(self.scdata.var.index)]
        self.genelist = self.genelist.sort_values()
        # Filter out genes not in gene list
        self.scdata = self.scdata[:,self.scdata.var.index.isin(self.genelist)]
        sc.pp.normalize_total(self.scdata, target_sum=1e6) # normalize to sum to 1,000,000
        # sc.pp.regress_out(scdata, ['total_counts'], n_jobs=1)
        # Transpose, filter out genes not in gene list, then sort column (by gene name)
        self.bkdata = self.bkdata.T.loc[:,self.genelist].sort_index(axis=1)

    def load_scdata(data_directories, cell_types):
        # Read and merge 10X Genomics scRNA-seq data
        scdata = None
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
        sns.jointplot(x='total_counts', y='n_genes_by_counts', height=8, data=scdata.obs,
            kind='scatter', hue='celltype')
        sns.jointplot(x='total_counts', y='pct_counts_mito', height=8, data=scdata.obs,
            kind='scatter', hue='celltype')
        sns.jointplot(x='total_counts', y='pct_counts_ribo', height=8, data=scdata.obs,
            kind='scatter', hue='celltype')
        sns.jointplot(x='total_counts', y='pct_counts_mribo', height=8, data=scdata.obs,
            kind='scatter', hue='celltype')
        plt.show()
        # Filter out cells with >5% of counts from mitochondria and mitoribosome
        scdata = scdata[scdata.obs.pct_counts_mito < 5, :]
        scdata = scdata[scdata.obs.pct_counts_mribo < 1, :]
        return scdata

    def __call__(self, whichdata):
        if whichdata == 'scdata':
            out = []
            for c in tqdm(self.celltypes):
                out.append(self.scdata[self.scdata.obs.celltype==c].to_df().sort_index(axis=1))
        elif whichdata == 'bkdata':
            out = self.bkdata
        elif whichdata == 'genelist':
            out = self.genelist
        else:
            raise ValueError('Choose only one of the following: "scdata", "bkdata", or "genelist"')
        return out
