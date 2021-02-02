import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns

data_directories = ["DATA/b_cells_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "DATA/cd4_t_helper_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "DATA/cd14_monocytes_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "DATA/cd34_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "DATA/cd56_nk_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "DATA/cytotoxic_t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "DATA/memory_t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "DATA/naive_cytotoxic_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "DATA/naive_t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                    "DATA/regulatory_t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/"]
cell_types = ['B_cell','CD4_helper','CD14','CD34','CD56_NK','CD8_cytotoxic','CD4_CD45RO_memory','CD8_CD45RA_naive','CD4_CD45RA_naive','CD4_CD25_regulatory']
bkdata_path = 'DATA/TCGA/TCGA_GDC_HTSeq_Counts.txt'
gene_list_path = 'DATA/Immune Gene Lists/genes.csv'

class DataPreprocess():
    def __init__(self, datadir, celltypes, bkdata_path, gene_list_path):
        self.datadir = datadir
        self.celltypes = celltypes
        self.scdata = load_scdata(self.datadir, self.celltypes)
        self.bkdata = pd.read_csv(bkdata_path, sep='\t')
        if gene_list_path is not None:
            self.genelist = pd.read_csv(gene_list_path, squeeze=True, names=['name'])
        else:
            self.genelist = None

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

    def filter_genes(bkdata, scdata, genelist):
        if genelist is not None:
            gene_list = genelist.drop_duplicates()
            gene_list = gene_list[gene_list.isin(bkdata.index.drop_duplicates())]
        else:
            gene_list = bkdata.index.drop_duplicates()
        gene_list = gene_list.sort_values()
        scdata = scdata[:,scdata.var.index.isin(gene_list)]
        return scdata

    def __call__(self):
        raw = self.filter_genes(self.bkdata, self.scdata, self.genelist)
        sc.pp.normalize_total(raw, target_sum=1e6) # normalize to sum to 1,000,000
        # sc.pp.regress_out(scdata, ['total_counts'], n_jobs=1)
        scdata = []
        for c in tqdm(self.celltypes):
            scdata.append(raw[raw.obs.celltype==c].to_df().sort_index(axis=1))
        return scdata
