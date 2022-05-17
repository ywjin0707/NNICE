library(reshape2)
library(tidyverse)
library(ggpubr)
meltextract <- function(res, method, sim_labels){
    if(!identical(colnames(res), label_cell_types)){
        colnames(res) <- label_cell_types
    }
    tmp <- reshape2::melt(res)
    tmp$method <- method
    tmp$true <- reshape2::melt(sim_labels)$value
    return(tmp)
}
fun_breaks = function(limits) {
    breaks = signif(max(limits) * c(0.25, 0.5, 0.75),1)
    names(breaks) = attr(breaks, "labels")
    breaks
}
fun_limits = function(limits) {
    lim = signif(max(limits) * c(0, 1.5),1)
    names(lim) = attr(lim, "labels")
    lim
}
###### ABIS monaco results ######
library(immunedeconv)
set_cibersort_binary('../../Preliminary Work/Rcode/CIBERSORT.R')
set_cibersort_mat('./LM22.txt')

label_cell_types <- c("B cell", "Monocytic lineage", "NK cell", "T cell CD4+", "T cell CD4+ (non-regulatory)", "T cell CD4+ naive", 
                      "T cell CD4+ memory", "T cell regulatory (Tregs)", "T cell CD8+", "T cell CD8+ naive", "T cell CD8+ memory",
                      "other cell")
label_cell_types <- c("B cell", "T cell CD4+", "T cell CD8+", "NK cell", "myeloid", "other cell")

label_abis <- read.csv("../../DATA/GSE107011/abis_pbmc_proportions.csv", sep = ",", header = TRUE, check.names = FALSE)
label_abis <- label_abis/100
label_abis <- data.frame(label_abis) %>%
    transmute("B cell" = B_cells,
              "T cell CD4+" = T_CD4,
              "T cell CD8+" = T_CD8,
              "NK cell" = NK,
              "Monocytic lineage" = Monocytes,
              "other cell" = other)

res_abis_tile_sim <- read.csv("../keras/logs/202007151522_ind_reg_rmseloss_l2001_sim_ABIS_2997_pred.csv", sep = ",", header = FALSE)
res_abis_tile_sim <- t(apply(t(res_abis_tile_sim), 2, function(x){x/sum(x)}))
res_abis_tile_real <- read.csv("../keras/logs/202007151522_ind_reg_rmseloss_l2001_real_ABIS_2997_pred.csv", sep = ",", header = FALSE)
res_abis_tile_real <- t(apply(t(res_abis_tile_real), 2, function(x){x/sum(x)}))
res_abis_timer <- deconvolute(abis, method = "timer", indications = c(rep("brca", ncol(abis))))
res_abis_cibersort <- deconvolute(abis, method = "cibersort")
res_abis_cibersort_abs <- deconvolute(abis, method = "cibersort_abs")
res_abis_quantiseq <- deconvolute(abis, method = "quantiseq", tumor = TRUE)
res_abis_epic <- deconvolute(abis, method = "epic")
# res_abis_xcell <- deconvolute(abis, method = "xcell")
res_abis_mcpcounter <- deconvolute(abis, method = "mcp_counter")

# res_SDY67_tile_sim <- read.csv("../keras/logs/202007151522_ind_reg_rmseloss_l2001_sim_SDY67_2997_pred.csv", sep = "," , header = FALSE)
# res_SDY67_tile_sim <- t(apply(t(res_SDY67_tile_sim), 2, function(x){x/sum(x)}))
# res_SDY67_tile_real <- read.csv("../keras/logs/20200731/202007311236_SDY67_3884_all_predictions.csv", sep = "," , header = FALSE)
# res_SDY67_tile_real <- t(apply(t(res_SDY67_tile_real), 2, function(x){x/sum(x)}))
res_SDY67_tile_real <- read.csv("../keras/logs/20201201/SDY67_predictions.csv", sep = ",", header = TRUE)
res_SDY67_tile_real <- res_SDY67_tile_real %>% filter(
    (Celltype == 'B' & Q == 0.75) | (Celltype == 'T_CD4' & Q == 0.5) | (Celltype == 'T_CD8' & Q ==0.5) | (Celltype == 'NK' & Q == 0.5) |
        (Celltype == 'Monocyte' & Q == 0.5) | (Celltype == 'Other' & Q == 0.5)) %>% select(Celltype, Predicted.proportion, True.proportion)
res_SDY67_tile_real_pred <- data.frame("B" = res_SDY67_tile_real$Predicted.proportion[which(res_SDY67_tile_real$Celltype == "B")],
                                  "T_CD4" = res_SDY67_tile_real$Predicted.proportion[which(res_SDY67_tile_real$Celltype == "T_CD4")],
                                  "T_CD8" = res_SDY67_tile_real$Predicted.proportion[which(res_SDY67_tile_real$Celltype == "T_CD8")],
                                  "NK" = res_SDY67_tile_real$Predicted.proportion[which(res_SDY67_tile_real$Celltype == "NK")],
                                  "Monocyte" = res_SDY67_tile_real$Predicted.proportion[which(res_SDY67_tile_real$Celltype == "Monocyte")],
                                  "Other" = res_SDY67_tile_real$Predicted.proportion[which(res_SDY67_tile_real$Celltype == "Other")])
res_SDY67_tile_real_true <- data.frame("B" = res_SDY67_tile_real$True.proportion[which(res_SDY67_tile_real$Celltype == "B")],
                                  "T_CD4" = res_SDY67_tile_real$True.proportion[which(res_SDY67_tile_real$Celltype == "T_CD4")],
                                  "T_CD8" = res_SDY67_tile_real$True.proportion[which(res_SDY67_tile_real$Celltype == "T_CD8")],
                                  "NK" = res_SDY67_tile_real$True.proportion[which(res_SDY67_tile_real$Celltype == "NK")],
                                  "Monocyte" = res_SDY67_tile_real$True.proportion[which(res_SDY67_tile_real$Celltype == "Monocyte")],
                                  "Other" = res_SDY67_tile_real$True.proportion[which(res_SDY67_tile_real$Celltype == "Other")])
res_SDY67_tile_real_pred <- t(apply(t(res_SDY67_tile_real_pred), 2, function(x){x/sum(x)}))

res_SDY67_timer <- deconvolute(SDY67, method = "timer", indications = c(rep("brca", ncol(SDY67))))
res_SDY67_cibersort <- deconvolute(SDY67, method = "cibersort")
res_SDY67_cibersort_abs <- deconvolute(SDY67, method = "cibersort_abs")
res_SDY67_quantiseq <- deconvolute(SDY67, method = "quantiseq", tumor = TRUE)
res_SDY67_epic <- deconvolute(SDY67, method = "epic")
# res_SDY67_xcell <- deconvolute(SDY67, method = "xcell")
res_SDY67_mcpcounter <- deconvolute(SDY67, method = "mcp_counter")

bulk_samples_test <- read.csv('../../DATA/simbulk/20210105/20210105_N-100_C-500_rep-1_simbulk_data.csv', header = FALSE, row.names = immport)
Nprop_test <- read.csv('../../DATA/simbulk/20210105/20210105_N-100_C-500_rep-1_simbulk_label.csv', header = FALSE, col.names = label_cell_types)

# res_simbulk_tile_sim <- read.csv("../keras/logs/202007151522_ind_reg_rmseloss_l2001_sim_simbulk_test_pred.csv", sep = "," , header = FALSE)
res_simbulk_tile_sim_pred <- read.csv("../keras/logs/20210105/simulated_predictions.csv", header=TRUE)
res_simbulk_tile_sim_pred <- t(apply(t(res_simbulk_tile_sim_pred), 2, function(x){x/sum(x)}))
res_simbulk_tile_sim_true <- read.csv("../keras/logs/20210105/simulated_true.csv", header=TRUE)
res_simbulk_tile_sim_true <- t(apply(t(res_simbulk_tile_sim_true), 2, function(x){x/sum(x)}))
# res_simbulk_tile_sim <- t(apply(t(res_simbulk_tile_sim), 2, function(x){x/sum(x)}))
# res_simbulk_tile_real <- read.csv("../keras/logs/202007151522_ind_reg_rmseloss_l2001_real_simbulk_test_pred.csv", sep = "," , header = FALSE)
# res_simbulk_tile_real <- t(apply(t(res_simbulk_tile_real), 2, function(x){x/sum(x)}))
# res_simbulk_tile_real <- read.csv('../keras/logs/20201203/simbulk_predictions.csv', sep = ',', header = TRUE, row.names = 1)
# res_simbulk_tile_real <- t(apply(t(res_simbulk_tile_real), 2, function(x){x/sum(x)}))

res_simbulk_timer <- deconvolute(bulk_samples_test, method = "timer", indications = c(rep("brca", ncol(bulk_samples_test))))
res_simbulk_cibersort <- deconvolute(bulk_samples_test, method = "cibersort")
res_simbulk_cibersort_abs <- deconvolute(bulk_samples_test, method = "cibersort_abs")
res_simbulk_quantiseq <- deconvolute(bulk_samples_test, method = "quantiseq", tumor = TRUE)
res_simbulk_epic <- deconvolute(bulk_samples_test, method = "epic")
# res_simbulk_xcell <- deconvolute(bulk_samples_test, method = "xcell")
res_simbulk_mcpcounter <- deconvolute(bulk_samples_test, method = "mcp_counter")

# colnames(res_abis_tile) <- cell_types_list
# res_abis_tile <- data.frame(res_abis_tile) %>% 
#     transmute("B cell" = B_cell, 
#               "Monocytic lineage" = CD14_monocyte,
#               "NK cell" = CD56_NK,
#               "T cell CD4+" = CD4_T_helper + CD4_CD45RO_memory_T +CD4_CD45RA_naive_T, CD4_CD25_regulatory_T,
#               "T cell CD4+ (non-regulatory)" = CD4_T_helper + CD4_CD45RO_memory_T + CD4_CD45RA_naive_T, 
#               "T cell CD4+ naive" = CD4_CD45RA_naive_T,
#               "T cell CD4+ memory" = CD4_CD45RO_memory_T,
#               "T cell regulatory (Tregs)" = CD4_CD25_regulatory_T,
#               "T cell CD8+" = CD8_cytotoxic_T + CD8_CD45RA_naive_T,
#               "T cell CD8+ naive" = CD8_CD45RA_naive_T, 
#               "T cell CD8+ memory" = CD8_cytotoxic_T,
#               "other cell" = CD34,
#               CD4_CD25_regulatory_T = NULL)

change_colnames <- function(x){
    require(dplyr)
    res <- data.frame(x) %>%
        transmute("B cell" = B.cell,
                  "T cell CD4+" = T.CD4,
                  "T cell CD8+" = T.CD8,
                  "NK cell" = NK,
                  "Monocytic lineage" = Monocytes,
                  "other cell" = other)
    return(res)
}
SDY67_label <- change_colnames(SDY67_label)

colnames(res_abis_tile_sim) <- colnames(SDY67_label)
# res_abis_tile_sim <- change_colnames(res_abis_tile_sim)
colnames(res_abis_tile_real) <- colnames(SDY67_label)
# res_abis_tile_real <- change_colnames(res_abis_tile_real)

colnames(res_SDY67_tile_sim) <- colnames(SDY67_label)
# res_SDY67_tile_sim <- change_colnames(res_SDY67_tile_sim)
colnames(res_SDY67_tile_real_pred) <- colnames(SDY67_label)
colnames(res_SDY67_tile_real_true) <- colnames(SDY67_label)
# res_SDY67_tile_real <- change_colnames(res_SDY67_tile_real)


colnames(res_simbulk_tile_sim) <- colnames(SDY67_label)
# res_simbulk_tile_sim <- change_colnames(res_simbulk_tile_sim)
colnames(res_simbulk_tile_real) <- label_cell_types
# res_simbulk_tile_real <- change_colnames(res_simbulk_tile_real) 
colnames(res_simbulk_tile_sim_pred) <- label_cell_types
colnames(res_simbulk_tile_sim_true) <- label_cell_types

colnames(Nprop_test) <- colnames(SDY67_label)
# Nprop_test <- change_colnames(Nprop_test)



tmp1 <- rbind(meltextract(data.frame(res_abis_tile_sim, check.names = FALSE), "ANN-SIM", label_abis), 
             meltextract(data.frame(res_abis_tile_real, check.names = FALSE), "ANN-REAL", label_abis), 
             meltextract(data.frame(t(map_result_to_celltypes(res_abis_cibersort_abs, use_cell_types = label_cell_types)), check.names = FALSE), "CBA", label_abis),
             meltextract(data.frame(t(map_result_to_celltypes(res_abis_cibersort, use_cell_types = label_cell_types)), check.names = FALSE), "CBS", label_abis),
             meltextract(data.frame(t(map_result_to_celltypes(res_abis_epic, use_cell_types = label_cell_types)), check.names = FALSE), "EPC", label_abis),
             meltextract(data.frame(t(map_result_to_celltypes(res_abis_mcpcounter, use_cell_types = label_cell_types)), check.names = FALSE), "MCP", label_abis),
             meltextract(data.frame(t(map_result_to_celltypes(res_abis_quantiseq, use_cell_types = label_cell_types)), check.names = FALSE), "QTS", label_abis),
             meltextract(data.frame(t(map_result_to_celltypes(res_abis_timer, use_cell_types = label_cell_types)), check.names = FALSE), "TMR", label_abis))

tmp2 <- rbind(#meltextract(data.frame(res_SDY67_tile_sim, check.names = FALSE), "ANN-SIM", SDY67_label), 
             meltextract(data.frame(res_SDY67_tile_real_pred, check.names = FALSE), "DQR", res_SDY67_tile_real_true), 
             meltextract(data.frame(t(map_result_to_celltypes(res_SDY67_cibersort_abs, use_cell_types = label_cell_types)), check.names = FALSE), "CBA", SDY67_label),
             meltextract(data.frame(t(map_result_to_celltypes(res_SDY67_cibersort, use_cell_types = label_cell_types)), check.names = FALSE), "CBS", SDY67_label),
             meltextract(data.frame(t(map_result_to_celltypes(res_SDY67_epic, use_cell_types = label_cell_types)), check.names = FALSE), "EPC", SDY67_label),
             meltextract(data.frame(t(apply(data.frame(t(map_result_to_celltypes(res_SDY67_mcpcounter, use_cell_types = label_cell_types)), check.names = FALSE), 1, function(x){x/sum(x, na.rm = TRUE)}))), "MCP", SDY67_label),
             meltextract(data.frame(t(map_result_to_celltypes(res_SDY67_quantiseq, use_cell_types = label_cell_types)), check.names = FALSE), "QTS", SDY67_label),
             meltextract(data.frame(t(map_result_to_celltypes(res_SDY67_timer, use_cell_types = label_cell_types)), check.names = FALSE), "TMR", SDY67_label))

tmp3 <- rbind(meltextract(data.frame(res_simbulk_tile_sim_pred, check.names = FALSE), "DQR", res_simbulk_tile_sim_true),
             # meltextract(data.frame(res_simbulk_tile_real, check.names = FALSE), "DQR", Nprop_test), 
             meltextract(data.frame(t(map_result_to_celltypes(res_simbulk_cibersort_abs, use_cell_types = label_cell_types)), check.names = FALSE), "CBA", Nprop_test),
             meltextract(data.frame(t(map_result_to_celltypes(res_simbulk_cibersort, use_cell_types = label_cell_types)), check.names = FALSE), "CBS", Nprop_test),
             meltextract(data.frame(t(map_result_to_celltypes(res_simbulk_epic, use_cell_types = label_cell_types)), check.names = FALSE), "EPC", Nprop_test),
             meltextract(data.frame(t(apply(data.frame(t(map_result_to_celltypes(res_simbulk_mcpcounter, use_cell_types = label_cell_types)), check.names = FALSE), 1, function(x){x/sum(x, na.rm = TRUE)}))),"MCP", Nprop_test),
             meltextract(data.frame(t(map_result_to_celltypes(res_simbulk_quantiseq, use_cell_types = label_cell_types)), check.names = FALSE), "QTS", Nprop_test),
             meltextract(data.frame(t(map_result_to_celltypes(res_simbulk_timer, use_cell_types = label_cell_types)), check.names = FALSE), "TMR", Nprop_test))

tmp3$method <- factor(tmp3$method, levels = c("DQR", "QTS", "EPC", "TMR", "MCP", "CBA", "CBS"))

ggplot(data = tmp3, aes(x=true, y=value)) + 
    geom_point(size=.2, alpha=1, colour = "grey") +
    # geom_text() +
    stat_smooth(color="black", method = "lm", size = 0.4, fullrange=TRUE) +
    stat_cor(method = "pearson", size = 2.5, color = "black") +
    facet_grid(method~variable, scales = "free", margins = "variable") +
    # scale_color_manual() +
    # scale_x_continuous(breaks = c(0.25, 0.5, 0.75), limits = c(0.0, 1.0)) +
    scale_y_continuous(breaks = fun_breaks) +
    scale_y_continuous(breaks = c(0.25, 0.5, 0.75), limits = c(0.0, 1.0)) + 
    theme(legend.position = "none", strip.text=element_text(size=8), panel.spacing = unit(.5, "mm"), axis.text = element_text(size=8)) +
    ylab("estimated fraction") +
    xlab("true fraction") +
    theme_bw()

#############################################################

rmse_loss <- function(y_true, y_pred){
    return(sqrt(mean((y_true - y_pred)^2)))
}

SDY67_rmse_method = c()
for(m in unique(tmp2$method)){
    tmp = filter(tmp2, tmp2$method == m)
    append(SDY67_rmse_method, rmse_loss(tmp$true, tmp$value))
}

SDY67_rmse_method_celltype = data.frame(array(0, dim = c(length(unique(tmp2$variable)), length(unique(tmp2$method))), dimnames = list(unique(tmp2$variable), unique(tmp2$method))))
for(c in label_cell_types){''
    for(m in tmp2$method){
        tmp <- tmp2 %>% filter(tmp2$method == m, tmp2$variable == c)
        SDY67_rmse_method_celltype[c,m] = rmse_loss(tmp$true, tmp$value)
    }
}
SDY67_rmse_method_celltype
#############################################################

perf_res <- read.csv("../keras/20200731_performance_results_simbulk_model_param.csv")
perf_res <- perf_res[,-1]
perf_res$kernel_constraint <- factor(perf_res$kernel_constraint)
perf_res$bias_constraint <- factor(perf_res$bias_constraint)
perf_res$lambda <- factor(perf_res$lambda)

ggplot(data = perf_res, aes(x = lambda, y = X2_rmse)) +
    geom_boxplot() +
    facet_grid(cols = vars(kernel_constraint), rows = vars(bias_constraint))

################################### 20210615
library(immunedeconv)
set_cibersort_binary('../../Preliminary Work/Rcode/CIBERSORT.R')
set_cibersort_mat('./LM22.txt')

label_cell_types <- c("B cell", "T cell CD4+", "T cell CD8+", "NK cell", "Monocytic lineage", "other cell")

abis_x <- read.csv('../../DATA/GSE107011/abis_pbmc_data.csv', sep=',', header=TRUE, check.names=FALSE, row.names=1)
abis_y <- read.csv('../../DATA/GSE107011/abis_pbmc_label_6ct.csv',sep=',',header=TRUE, row.names = 1, check.names = FALSE)
sdy67_x <- read.csv('../../DATA/SDY67/SDY67_468_pp.csv', sep=',', header=TRUE,row.names=1, check.names=FALSE)
sdy67_y <- read.csv('../../DATa/SDY67/SDY67_468_label_6ct_pp.csv', sep=',', header=TRUE, row.names=1, check.names=FALSE)
simbulk_x <- read.csv('../../DATA/simbulk/20201105_N-1000_C-500_simbulk_data.csv', sep=',', header=FALSE)
simbulk_y <- read.csv('../../DATA/simbulk/20201105_N-1000_C-500_simbulk_label.csv', sep=',', header=FALSE)
rownames(simbulk_x) <- rownames(mydata2@assays$SCT)
colnames(simbulk_y) <- label_cell_types

### ABIS results
res_abis_nnice_sim <- read.csv("../../dice/log/predictions/20210614_230656_simbulk_ABIS_avg_predictions_6ct.csv", sep = ",", header = TRUE)
colnames(res_abis_nnice_sim) <- label_cell_types
res_abis_nnice_real <- read.csv("../../dice/log/predictions/20210614_235958_SDY67_ABIS_avg_predictions_6ct.csv", sep = ",", header = TRUE)
colnames(res_abis_nnice_real) <- label_cell_types
res_abis_timer <- deconvolute(abis_x, method = "timer", indications = c(rep("dlbc", ncol(abis_x))))
res_abis_cibersort <- deconvolute(abis_x, method = "cibersort")
# res_abis_cibersort_abs <- deconvolute(abis_x, method = "cibersort_abs")
res_abis_quantiseq <- deconvolute(abis_x, method = "quantiseq", tumor = FALSE)
res_abis_epic <- deconvolute(abis_x, method = "epic")
# res_abis_xcell <- deconvolute(abis, method = "xcell")
res_abis_mcpcounter <- deconvolute(abis_x, method = "mcp_counter")

tmp1 <- rbind(#meltextract(data.frame(res_abis_nnice_sim, check.names = FALSE), "NNICE-sim", abis_y), 
              meltextract(data.frame(res_abis_nnice_real, check.names = FALSE), "NNICE", abis_y), #, "NNICE-real", abis_y), 
              # meltextract(data.frame(t(apply(data.frame(t(map_result_to_celltypes(res_abis_cibersort_abs, use_cell_types = label_cell_types)), check.names = FALSE), 1, function(x){x/sum(x, na.rm = TRUE)}))), "CBA", abis_y),
              meltextract(data.frame(t(map_result_to_celltypes(res_abis_cibersort, use_cell_types = label_cell_types)), check.names = FALSE), "CBS", abis_y),
              meltextract(data.frame(t(map_result_to_celltypes(res_abis_epic, use_cell_types = label_cell_types)), check.names = FALSE), "EPC", abis_y),
              meltextract(data.frame(t(apply(data.frame(t(map_result_to_celltypes(res_abis_mcpcounter, use_cell_types = label_cell_types)), check.names = FALSE), 1, function(x){x/sum(x, na.rm = TRUE)}))), "MCP", abis_y),
              meltextract(data.frame(t(map_result_to_celltypes(res_abis_quantiseq, use_cell_types = label_cell_types)), check.names = FALSE), "QTS", abis_y),
              meltextract(data.frame(t(map_result_to_celltypes(res_abis_timer, use_cell_types = label_cell_types)), check.names = FALSE), "TMR", abis_y))

# tmp1$method <- factor(tmp1$method, levels = c("NNICE-sim", "NNICE-real", "QTS", "CBA", "EPC",  "MCP", "TMR",  "CBS"))
tmp1$method <- factor(tmp1$method, levels = c("NNICE", "QTS", "CBA", "EPC",  "MCP", "TMR",  "CBS"))

ggplot(data = tmp1, aes(x=true, y=value)) + 
    geom_point(size=0.5, colour = "grey") +
    # geom_text() +
    stat_smooth(color="black", method = "lm", size = 0.4, fullrange = TRUE) +
    stat_cor(method = "pearson", size = 3.5, color = "black", fontface='bold') +
    facet_grid(method~variable, scales = "free", margins = "variable") +
    # scale_color_manual() +
    scale_x_continuous(breaks = c(0.25, 0.5, 0.75), limits = c(0.0, 1.0)) +
    scale_y_continuous(breaks = fun_breaks) +
    scale_y_continuous(breaks = c(0.25, 0.5, 0.75), limits = c(0.0, 1.0)) + 
    theme(legend.position = "none", strip.text=element_text(size=8), panel.spacing = unit(.5, "mm"), axis.text = element_text(size=8)) +
    ylab("estimated fraction") +
    xlab("true fraction") +
    theme_bw()
ggsave(filename ='C:/Users/yw_ji/Documents/MSc Thesis/WRITING/perf_compare_ABIS_noNNICEsim2.png', width=12, height=9, device='png', dpi=600)

### SDY67 results
res_SDY67_nnice_sim <- read.csv("../../dice/log/predictions/20210614_230656_simbulk_SDY67_avg_predictions_6ct.csv", sep = ",", header = TRUE)
colnames(res_SDY67_nnice_sim) <- label_cell_types
res_SDY67_nnice_real <- read.csv("../../dice/log/predictions/20210614_210023_SDY67_10CV_avg_predictions_6ct.csv", sep = ",", header = TRUE)
colnames(res_SDY67_nnice_real) <- label_cell_types
res_SDY67_timer <- deconvolute(sdy67_x, method = "timer", indications = c(rep("dlbc", ncol(sdy67_x))))
res_SDY67_cibersort <- deconvolute(sdy67_x, method = "cibersort")
# res_SDY67_cibersort_abs <- deconvolute(sdy67_x, method = "cibersort_abs")
res_SDY67_quantiseq <- deconvolute(sdy67_x, method = "quantiseq", tumor = FALSE)
res_SDY67_epic <- deconvolute(sdy67_x, method = "epic")
# res_SDY67_xcell <- deconvolute(sdy67_x, method = "xcell")
res_SDY67_mcpcounter <- deconvolute(sdy67_x, method = "mcp_counter")

tmp2 <- rbind(#meltextract(data.frame(res_SDY67_nnice_sim, check.names = FALSE), "NNICE-sim", sdy67_y), 
              meltextract(data.frame(res_SDY67_nnice_real, check.names = FALSE), "NNICE", sdy67_y),#"NNICE-real", sdy67_y), 
              # meltextract(data.frame(t(map_result_to_celltypes(res_SDY67_cibersort_abs, use_cell_types = label_cell_types)), check.names = FALSE), "CBA", sdy67_y),
              meltextract(data.frame(t(map_result_to_celltypes(res_SDY67_cibersort, use_cell_types = label_cell_types)), check.names = FALSE), "CBS", sdy67_y),
              meltextract(data.frame(t(map_result_to_celltypes(res_SDY67_epic, use_cell_types = label_cell_types)), check.names = FALSE), "EPC", sdy67_y),
              meltextract(data.frame(t(apply(data.frame(t(map_result_to_celltypes(res_SDY67_mcpcounter, use_cell_types = label_cell_types)), check.names = FALSE), 1, function(x){x/sum(x, na.rm = TRUE)}))), "MCP", sdy67_y),
              meltextract(data.frame(t(map_result_to_celltypes(res_SDY67_quantiseq, use_cell_types = label_cell_types)), check.names = FALSE), "QTS", sdy67_y),
              meltextract(data.frame(t(map_result_to_celltypes(res_SDY67_timer, use_cell_types = label_cell_types)), check.names = FALSE), "TMR", sdy67_y))

# tmp2$method <- factor(tmp2$method, levels = c("NNICE-sim", "NNICE-real", "QTS", "CBA", "EPC",  "MCP", "TMR",  "CBS"))
tmp2$method <- factor(tmp2$method, levels = c("NNICE", "QTS", "CBA", "EPC",  "MCP", "TMR",  "CBS"))


ggplot(data = tmp2, aes(x=true, y=value)) + 
    geom_point(size=.5, alpha=1, colour = "grey") +
    # geom_text() +
    stat_smooth(color="black", method = "lm", size = 0.4, fullrange = TRUE) +
    stat_cor(method = "pearson", size = 3.5, color = "black", fontface='bold') +
    facet_grid(method~variable, scales = "free", margins = "variable") +
    # scale_color_manual() +
    # scale_x_continuous(breaks = c(0.25, 0.5, 0.75), limits = c(0.0, 1.0)) +
    scale_y_continuous(breaks = fun_breaks) +
    scale_y_continuous(breaks = c(0.25, 0.5, 0.75), limits = c(0.0, 1.0)) + 
    theme(legend.position = "none", strip.text=element_text(size=8), panel.spacing = unit(.5, "mm"), axis.text = element_text(size=8)) +
    ylab("estimated fraction") +
    xlab("true fraction") +
    theme_bw()
ggsave(filename ='C:/Users/yw_ji/Documents/MSc Thesis/WRITING/perf_compare_SDY67_noNNICEsim2.png', width=12, height=9, device='png', dpi=600)

### simbulk results
res_simbulk_nnice_sim <- read.csv("../../dice/log/predictions/20210614_230656_simbulk_test_avg_predictions_6ct.csv", sep = ",", header = TRUE)
colnames(res_simbulk_nnice_sim) <- label_cell_types
res_simbulk_nnice_real <- read.csv("../../dice/log/predictions/20210614_235958_SDY67_simbulk_avg_predictions_6ct.csv", sep = ",", header = TRUE)
colnames(res_simbulk_nnice_real) <- label_cell_types
res_simbulk_timer <- deconvolute(simbulk_x, method = "timer", indications = c(rep("dlbc", ncol(simbulk_x))))
res_simbulk_cibersort <- deconvolute(simbulk_x, method = "cibersort")
# res_simbulk_cibersort_abs <- deconvolute(simbulk_x, method = "cibersort_abs")
res_simbulk_quantiseq <- deconvolute(simbulk_x, method = "quantiseq", tumor = FALSE)
res_simbulk_epic <- deconvolute(simbulk_x, method = "epic")
# res_simbulk_xcell <- deconvolute(simbulk_x, method = "xcell")
res_simbulk_mcpcounter <- deconvolute(simbulk_x, method = "mcp_counter")

tmp3 <- rbind(meltextract(data.frame(res_simbulk_nnice_sim, check.names = FALSE), "NNICE", simbulk_y),#, "NNICE-sim", simbulk_y), #res_simbulk_tile_sim_true
              # meltextract(data.frame(res_simbulk_nnice_real, check.names = FALSE), "NNICE-real", simbulk_y), # Nprop_test
              # meltextract(data.frame(t(map_result_to_celltypes(res_simbulk_cibersort_abs, use_cell_types = label_cell_types)), check.names = FALSE), "CBA", simbulk_y),
              meltextract(data.frame(t(map_result_to_celltypes(res_simbulk_cibersort, use_cell_types = label_cell_types)), check.names = FALSE), "CBS", simbulk_y),
              meltextract(data.frame(t(map_result_to_celltypes(res_simbulk_epic, use_cell_types = label_cell_types)), check.names = FALSE), "EPC", simbulk_y),
              meltextract(data.frame(t(apply(data.frame(t(map_result_to_celltypes(res_simbulk_mcpcounter, use_cell_types = label_cell_types)), check.names = FALSE), 1, function(x){x/sum(x, na.rm = TRUE)}))),"MCP", simbulk_y),
              meltextract(data.frame(t(map_result_to_celltypes(res_simbulk_quantiseq, use_cell_types = label_cell_types)), check.names = FALSE), "QTS", simbulk_y),
              meltextract(data.frame(t(map_result_to_celltypes(res_simbulk_timer, use_cell_types = label_cell_types)), check.names = FALSE), "TMR", simbulk_y))

# tmp3$method <- factor(tmp3$method, levels = c("NNICE-sim", "NNICE-real", "QTS", "CBA", "EPC",  "MCP", "TMR",  "CBS"))
tmp3$method <- factor(tmp3$method, levels = c("NNICE", "QTS", "CBA", "EPC",  "MCP", "TMR",  "CBS"))


ggplot(data = tmp3, aes(x=true, y=value)) + 
    geom_point(size=.5, alpha=1, colour = "grey") +
    # geom_text() +
    stat_smooth(color="black", method = "lm", size = 0.4, fullrange = TRUE) +
    stat_cor(method = "pearson", size = 3.5, color = "black", fontface='bold') +
    facet_grid(method~variable, scales = "free", margins = "variable") +
    # scale_color_manual() +
    # scale_x_continuous(breaks = c(0.25, 0.5, 0.75), limits = c(0.0, 1.0)) +
    scale_y_continuous(breaks = fun_breaks) +
    scale_y_continuous(breaks = c(0.25, 0.5, 0.75), limits = c(0.0, 1.0)) + 
    theme(legend.position = "none", strip.text=element_text(size=8), panel.spacing = unit(.5, "mm"), axis.text = element_text(size=8)) +
    ylab("estimated fraction") +
    xlab("true fraction") +
    theme_bw()
ggsave(filename ='C:/Users/yw_ji/Documents/MSc Thesis/WRITING/perf_compare_simbulk_noNNICEreal2.png', width=12, height=9, device='png', dpi=600)
