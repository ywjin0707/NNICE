library(Seurat)
library(patchwork)
library(ggplot2)

input_directories <- c("../../DATA/b_cells_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                       "../../DATA/cd4_t_helper_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                       "../../DATA/cd14_monocytes_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                       "../../DATA/cd34_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                       "../../DATA/cd56_nk_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                       "../../DATA/cytotoxic_t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                       "../../DATA/memory_t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                       "../../DATA/naive_cytotoxic_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                       "../../DATA/naive_t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/",
                       "../../DATA/regulatory_t_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/")

## Check if all scRNA-seq datasets have same rownames
# for(K in 2:length(input_directories)){print(identical(rownames(raw), rownames(Read10X(input_directories[K], gene.column =1))))}
## Yes, they do. All 10 datasets return """TRUE"""

cell_types_list <- c("B_cell", "CD4_T_helper", "CD14_monocyte", "CD34", "CD56_NK", "CD8_cytotoxic_T", "CD4_CD45RO_memory_T", "CD8_CD45RA_naive_T",
                     "CD4_CD45RA_naive_T", "CD4_CD25_regulatory_T")

### Load validation datasets ###
tcga_tpm <- read.csv("../../DATA/TCGA/TCGA_GDC_HTSeq_TPM.txt", sep = "\t")
load("C:/Users/yw_ji/Documents/MSc Thesis/DATA/EPIC/melanoma_data.rda") # This loads and object named: "melanoma_data"
abis <- read.csv("../../DATA/GSE107011/GSE107011_Processed_data_TPM.txt",sep = "\t")
abis$X <- gsub("\\..*", "", abis$X) # Ensembl IDs for this dataset has "." followed by version number, which will be removed


raw <- Read10X(data.dir = input_directories[1], gene.column = 2)
# Retrieve ensembl and hgnc ids for conversion later
# Here, we confirmed that all 10 datasets have identical rownames above so we don't need to worry about matching
hgnc_id <- rownames(raw) # retrieve hgnc_ids to match later
mit_id <- grep("^MT-", rownames(raw)) # extract positions of mitochondrial genes

# Initiate the merged Seurat object
raw <- Read10X(data.dir = input_directories[1], gene.column = 1) # using ensembl ids because it is more stable
ensembl_id <- rownames(raw) # retrieve ensembl_ids to match later
key_table <- data.frame(hgnc_id, ensembl_id, stringsAsFactors = FALSE)
# Add cell type identity to cell barcode IDs (colnames)
colnames(raw) <- paste(cell_types_list[1], colnames(raw), sep = "_")

# Not all variable genes selected by SCTransform will be found in validation datasets. 
# Therefore, only those genes/features that are in validation datasets will be considered. 
library(AnnotationDbi)
library(EnsDb.Hsapiens.v79)
# need to convert hgnc id to ensembl id
genes_melanoma <- mapIds(EnsDb.Hsapiens.v79, keys=rownames(melanoma_data$counts), column = "GENEID", keytype = "SYMBOL") 
genes_tcga <- mapIds(EnsDb.Hsapiens.v79, keys=rownames(tcga_tpm), column = "GENEID", keytype = "SYMBOL")
genes_abis <- abis$X

genes_all <- rownames(raw)[rownames(raw) %in% abis$X]
genes_all <- genes_all[genes_all %in% genes_melanoma]
genes_all <- genes_all[genes_all %in% genes_tcga]
genes_all <- c(rownames(raw)[mit_id], genes_all) # add back in mitochondrial genes because they were filtered out
length(genes_all) # 20011 genes/features overlap across scRNAseq and bulk validation datasets + 13 mitochondrial genes = 20024 genes/features

convert_id <- match(genes_all, key_table$ensembl_id)
genes_all_hgnc <- key_table$hgnc_id[convert_id]
# Proceeding with this produces """Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')"""
# Find which HGNC symbols contain "_"
genes_all_hgnc[grep("_", genes_all_hgnc)] #8969 - "RP11-453F18__B.1"
genes_all_hgnc[8969] <- "RP11-453F18"

# Subseting genes/features must occur before creation of Seurat object
raw <- raw[genes_all,]
rownames(raw) <- genes_all_hgnc
obj <- CreateSeuratObject(raw, min.cells = 0, project = cell_types_list[1], 
                          min.features = 200, names.field = NULL, names.delim = NULL)

# Read in the rest of the datasets and merge
for(K in 2:length(input_directories)){
    raw <- Read10X(data.dir = input_directories[K], gene.column = 1) # using ensembl ids because it is more stable
    colnames(raw) <- paste(cell_types_list[K], colnames(raw), sep = "_")
    # raw <- raw[genes_all,]
    # rownames(raw) <- genes_all_hgnc
    raw <- CreateSeuratObject(raw, min.cells = 0, project = cell_types_list[K], 
                              min.features = 200, names.field = NULL, names.delim = NULL)
    obj <- merge(x = obj, y = raw) # merge all data sets one by one
}
obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt") # calculate percentage of mitochondrial genes
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
### ADD LATER CONSIDERATIONS FOR RIBOSOMAL GENES (pattern = "^RP[LS]") ###
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# normalize using negative binomial regression with regularization, accounting for mitochondrial genes.Slow - takes about 1 hour. #2:23 -
obj <- SCTransform(obj, vars.to.regress = "percent.mt", variable.features.n = 3000, conserve.memory = TRUE) 
obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt_SCT", assay = "SCT")

## Draw up a quality control plot to see results of QC
VlnPlot(obj, features = c("nCount_SCT", "nFeature_SCT", "percent.mt_SCT"), ncol = 3)
rm(raw,convert_id,ensembl_id, genes_abis, genes_tcga, genes_melanoma, hgnc_id, mit_id)

saveRDS(obj, file = "2020_06_08_obj_SCT_min0cell_min200max3000features_max5percentmt_3000varfeatures_conservemem.rds")

grep("^MT-", obj@assays$SCT@var.features) # 1537, 1944, 2810
mydata <- subset(obj, features = obj@assays$SCT@var.features[c(-1537, -1944, -2810)])

# output <- as.matrix(obj@assays$SCT[obj@assays$SCT@var.features]) # Keep only the 3000 variable genes
# id <- match(rownames(output), gsub("_", "-", gene_names))
# rownames(output) <- gene_names[id]
# colnames(output) <- gsub("_.{14}$", "", colnames(output))



### ABIS Monaco dataset of PBMC samples
# abis <- read.csv("../../DATA/GSE107011/GSE107011_Processed_data_TPM.txt", sep = "\t")
# abis$X <- gsub("\\..*", "", abis$X)
library(stringr)
abis <- abis[,c(1,which(str_sub(colnames(abis), -4, -1) == "PBMC"))] # only keep 12 PBMC samples
colnames(abis) <- substr(gsub("X", "", colnames(abis)), 1, 4)
abis <- abis[,-8] # exclude one without proportions data

id <- match(rownames(mydata), key_table$hgnc_id)
id <- match(key_table$ensembl_id[id], abis[,1])
abis_dnn <- abis
rownames(abis_dnn) <- abis_dnn[,1]
abis_dnn <- abis_dnn[,-1]
abis_dnn <- abis_dnn[id,]
write.csv(abis_dnn, "../../DATA/GSE107011/GSE107011_Processed_data_TPM_HGNC_Gene_ID_filtered_2997.csv", row.names = FALSE)

### use all of ABIS samples
abis <- read.csv("../../DATA/GSE107011/GSE107011_Processed_data_TPM.txt",sep = "\t", check.names = FALSE)
colnames(abis)[1] <- 'ID'
abis$ID <- gsub("\\..*", "", abis$ID) # Ensembl IDs for this dataset has "." followed by version number, which will be removed
library(AnnotationDbi)
library(EnsDb.Hsapiens.v79)
id <- mapIds(EnsDb.Hsapiens.v79, keys=abis$ID, column = "SYMBOL", keytype = "GENEID") 
abis$ID <- id
abis <- abis[complete.cases(abis),]
abis <- aggregate(.~ID, data = abis, max)
rownames(abis) <- abis$ID
abis <- abis[,-1]
write.csv(abis, "../../DATA/GSE107011/GSE107011_Processed_data_TPM_HGNC_Gene_ID.csv")

### 
library(AnnotationDbi)
library(EnsDb.Hsapiens.v79)
# need to convert ensembl id to hgnc id
id <- mapIds(EnsDb.Hsapiens.v79, keys=abis[,1], column = "SYMBOL", keytype = "GENEID") 
abis[,1] <- id
abis <- abis[complete.cases(abis),]
colnames(abis)[1] <- "hgnc"
abis <- aggregate(.~hgnc, data = abis, max)
rownames(abis) <- abis$hgnc
abis <- abis[,-1]

abis_immport <- abis
abis_immport <- abis_immport[rownames(abis_immport) %in% key_table$hgnc_id[match(immport, key_table$ensembl_id)],]
abis_immport <- abis_immport[rownames(abis_immport) %in% rownames(SDY67_immport),]
tmp <- data.frame(matrix(0, nrow = 5, ncol = ncol(abis_immport)), 
                  row.names = rownames(SDY67_immport)[which(!(rownames(SDY67_immport) %in% rownames(abis_immport)))])
colnames(tmp) <- colnames(abis_immport)
abis_immport <- rbind(abis_immport, tmp)
abis_immport <- abis_immport[match(rownames(SDY67_immport), rownames(abis_immport)),]
write.csv(abis_immport, "../../DATA/GSE107011/GSE107011_Processed_data_TPM_HGNC_Gene_ID_filtered_3884.csv", row.names = FALSE)

####### SDY67 Zimmerman ########

SDY67_1 <- read.csv("../../DATA/SDY67/SDY67_EXP14625_RNA_seq.703317.tsv", sep = "\t")
SDY67_2 <- read.csv("../../DATA/SDY67/SDY67_EXP13377_RNA_seq.703318.tsv", sep = "\t")
id <- match(SDY67_1$GENE_SYMBOL, SDY67_2$GENE_SYMBOL)
SDY67 <- cbind(SDY67_1, SDY67_2[id,-1])
rm(SDY67_1, SDY67_2)
rownames(SDY67) <- SDY67$GENE_SYMBOL
SDY67 <- SDY67[,-1]

SDY67_id <- read.table("../../DATA/SDY67/SDY67-DR34_Subject_2_RNA_sequencing_result.txt", sep = "\t", header = TRUE)
id <- match(colnames(SDY67), SDY67_id$Expsample.Accession)
colnames(SDY67) <- paste(SDY67_id[id,]$Subject.Accession, SDY67_id[id,]$Study.Time.Collected, sep = "_")

SDY67_label <- read.csv("../../DATA/SDY67/SDY67_extracted_from_mmc7.csv", sep = ",", header = TRUE, row.names = 1, check.names = FALSE)
SDY67_label <- SDY67_label[complete.cases(SDY67_label),]
SDY67_label$other <- 0
SDY67_label <- SDY67_label/100
# SDY67_label[is.na(SDY67_label)] <- 0
SDY67_label <- apply(SDY67_label, 1, function(x){
    if(sum(x) >= 1){
        x <- x/sum(x)
    } else {
        x[6] <- (1-sum(x))
        x <- x/sum(x)
    }
})

SDY67_label <- t(SDY67_label)

SDY67 <- SDY67[,colnames(SDY67) %in% rownames(SDY67_label)]
SDY67_dnn <- SDY67
id <- match(rownames(mydata), rownames(SDY67_dnn))
SDY67_dnn <- SDY67_dnn[id,]
rownames(SDY67_dnn) <- rownames(mydata)
SDY67_dnn[is.na(SDY67_dnn)] <- 0
write.csv(SDY67_dnn, "../../DATA/SDY67/SDY67_2997.csv", row.names = FALSE)


id <- match(colnames(SDY67), rownames(SDY67_label))
SDY67_label <- SDY67_label[id,]
SDY67_label <- data.frame(SDY67_label[complete.cases(SDY67_label),])
write.csv(SDY67_label, "../../DATA/SDY67/SDY67_label.csv", row.names = FALSE)

####### STAT-CAN BREAST CANCER DATA #######

incidence <- c(0,0,5,0,0,5,5,
             100,90,85,100,95,90,105,
             705,705,635,725,780,750,795,
             2740,2625,2545,2585,2565,2545,2590,
             4220,4155,4205,4320,4485,4580,4410,
             4585,4555,4855,5100,5070,5245,5185,
             3175,3225,3385,3360,3525,3745,3805,
             1895,1805,1780,1825,1795,1870,1875,
             385,405,445,460,450,515,400)

incidence <- matrix(incidence, ncol = 7, byrow = TRUE)
rownames(incidence) <- c("0-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+")
colnames(incidence) <- as.character(2011:2017)

incidence_prop <- apply(incidence, 2, function(x){x/sum(x)})

library(ggplot2)
library(reshape2)
tmp <- melt(incidence_prop)
colnames(tmp) <- c("Age Group", "Year", "Proportion")
tmp$Year <- as.factor(tmp$Year)
g <- ggplot(data = tmp, aes(Year, Proportion)) + 
    geom_bar(stat = "identity", aes(fill=`Age Group`), position=position_fill(reverse=TRUE)) + scale_y_continuous(breaks = seq(0,1, by=0.1)) +
    scale_fill_brewer(palette="Spectral") + theme_bw()
g

death <- read.csv("../../DATA/CCR/females_total_deaths_mod.csv", row.names = 1, header=TRUE)
colnames(death) <- gsub("X", "", colnames(death))


#### Immport gene list ####

immport <- read.csv("../../DATA/Immune Gene Lists/immport.csv")
immport <- as.character(immport$ensembl)
immport <- immport[immport %in% key_table$ensembl_id]
SDY67_immport <- SDY67
SDY67_immport <- SDY67_immport[rownames(SDY67_immport) %in% key_table$hgnc_id[match(immport, key_table$ensembl_id)],]
write.csv(SDY67_immport, "../../DATA/SDY67/SDY67_3884.csv", row.names = FALSE)

### Filter genes from TCGA & METABRIC ###
tcga_tpm_2925 <- tcga_tpm[match(rownames(mydata2),rownames(tcga_tpm)),]
write.csv(tcga_tpm_2925, '../../DATA/TCGA/TCGA_TPM_2925.csv', row.names = FALSE)

metabric <- read.csv('../../DATA/METABRIC/data_expression_median.txt', sep='\t')
rownames(metabric) <- metabric$Hugo_Symbol
metabric <- metabric[,c(-1,-2)]
id <- match(rownames(mydata2),rownames(metabric))
metabric_2925 <- metabric[id,]
metabric_2925[is.na(metabric_2925)] <- 0 # Assign zero to all elements with NA rownames
rownames(metabric_2925)[which(is.na(id))] <- rownames(mydata2)[which(is.na(id))] # Assign rowname to unmatched rows
write.csv(metabric_2925, '../../DATA/METABRIC/METABRIC_2925.csv', row.names = FALSE)
# table(rownames(metabric_2925)=='NA') # There is only one rowname missing in METABRIC
# which(rownames(metabric_2925)=='NA')
# rownames(metabric_2925)[4] <- 'CDK11B'
# metabric_2925[4,] <- 0

summary(apply(bulk_samples_test, 2, function(x){summary(x)[4]}))
summary(apply(tcga_tpm_2925, 2, function(x){summary(x)[4]}))
summary(apply(metabric_2925, 2, function(x){summary(x)[4]}))
summary(apply(tmp, 2, function(x){summary(x)[4]}))
