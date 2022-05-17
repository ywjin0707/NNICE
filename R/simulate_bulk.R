library(MCMCpack)
library(dplyr)


make_bulk_sample <- function(seuratobj, celltype_ids, celltype_counts){
  cell_ids = lapply(names(celltype_ids), function(celltype){
    cell_count = celltype_counts[celltype]
    sample(celltype_ids[[celltype]], size=cell_count, replace=TRUE)
  })
  bulk_sample = rowSums(seuratobj[,unlist(cell_ids)])
  return(bulk_sample)
}



make_bulk_dataset <- function(seuratobj, celltypes_list, N=1000, C=500){
  celltype_ids <- list() # Initialize list
  for(K in 1:length(celltypes_list)){
    celltype_ids[[K]] <- WhichCells(seuratobj, ident = celltypes_list[[K]])
  }
  names(celltype_ids) <- names(celltypes_list)
  celltype_fractions <- rdirichlet(N, c(rep(1, times = length(simbulk_cell_types))))
  colnames(celltype_fractions) = names(celltypes_list)
  celltype_counts <- celltype_fractions * C
  expr = pblapply(1:nrow(celltype_counts), function(i){
    make_bulk_sample(seuratobj, celltype_ids, celltype_counts[i,])
  })
  expr = do.call(cbind, expr)
  colnames(expr) = 1:ncol(expr)
  rownames(expr) = rownames(seuratobj)
  rownames(celltype_fractions) = colnames(expr)
  
  return(list(expr, celltype_fractions))
  
}

# Example usage
# celltypes_list <- list("B_cell", c("CD4_CD45RA_naive_T", "CD4_CD45RO_memory_T", "CD4_T_helper", "CD4_CD25_regulatory_T"), 
#                        c("CD8_CD45RA_naive_T", "CD8_cytotoxic_T"), "CD56_NK", "CD14_monocyte", "CD34" )
# names(celltypes_list) <- c('B_cell','CD4_cell','CD8_cell','NK_cell','monocyte','other')
# mydata <- make_bulk_dataset(obj, simbulk_cell_types, N=10)
