from datetime import datetime
from model import *
from data import *


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
cell_types = ['B_cell','CD4_helper','CD14','CD34','CD56_NK','CD8_cytotoxic','CD4_CD45RO_memory','CD8_CD45RA_naive','CD4_CD45RA_naive','CD4_CD25_regulatory']
bkdata_path = '../DATA/TCGA/TCGA_GDC_HTSeq_TPM.csv'
# gene_list_path = '../DATA/Immune Gene Lists/genes.csv'
data_paths = ['../DATA/TCGA/TCGA_GDC_HTSeq_TPM.csv',
              '../DATA/METABRIC/METABRIC.csv',
              '../DATA/SDY67/SDY67_477.csv',
              '../DATA/Gene Lists/immport_genelist.csv',
              '../DATA/Gene Lists/scdata_genelist_filtered.csv']


if __name__ == '__main__':
    EPOCHS=1000
    features = FeatureList(data_paths)
    myData = DataPreprocess(data_directories, cell_types, bkdata_path, features)
    simulator = Build_Simulator(myData('scdata'))
    discriminator = Build_Discriminator()

    checkpoint_dir = './log/'
    checkpoint_prefix = os.path.join(checkpoint_dir, "checkpoint")
    checkpoint = tf.train.Checkpoint(Optmzr_sim=Optmzr_sim,
                                     Optmzr_dis=Optmzr_dis,
                                     simulator=simulator,
                                     discriminator=discriminator)

    myModel = train(myData('bkdata'), epochs=EPOCHS)
