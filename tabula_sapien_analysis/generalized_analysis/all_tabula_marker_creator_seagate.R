### The purpose of this file is to create markers for all the genes of interest

#####
### Load packages ###
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(dplyr)
library(stringr)
source("global_in/general_scripts.R")
#####
### Set up environment ###
tabula_key <- read.csv("global_in/general_tissues_tabula_mkh.csv")
gene_location <- "MARCC output/variance_genes_kendall/"
data_location <- "D:/tabula_sapiens_data/"
output_locations <- "tabula_sapien_analysis/generalized_analysis/data_out/"
#####
### Nested for loop on all tabula types
tabula_key <- filter(tabula_key, tabula_sapiens != "SKIP")
unique_tabula_tissues <- unique(tabula_key$tabula_sapiens)
## Remove "&" samples and clean up looping vector
unique_tabula_tissues <- unique_tabula_tissues[!grepl("&", unique_tabula_tissues)]
unique_tabula_tissues <- gsub(" ", "_",unique_tabula_tissues)
unique_tabula_tissues[unique_tabula_tissues == "All_tabula"] <- "TabulaSapiens"
unique_tabula_tissues <- unique_tabula_tissues[unique_tabula_tissues != "Small_intestine"]
## put TabulaSapiens last
unique_tabula_tissues <- c(unique_tabula_tissues[-which(unique_tabula_tissues == "TabulaSapiens")], "TabulaSapiens")

master_tables <- NULL

for (tabula_ind in seq_along(unique_tabula_tissues)) {
    cur_tabula <- unique_tabula_tissues[tabula_ind]
    if (cur_tabula == "TabulaSapiens") {
        file_of_interest <- paste0(data_location,"/",cur_tabula)  
    } else{
        file_of_interest <- paste0(data_location,"/TS_",cur_tabula)
    }
    ## Prepare seurat data
    SeuratDisk::Convert(paste0(file_of_interest,".h5ad"), dest = "h5seurat", overwrite = TRUE,assay = "RNA")
    cur_dat <- LoadH5Seurat(paste0(file_of_interest,".h5seurat"), assays = "RNA")
    ## Load tissue data
    if (cur_tabula == "TabulaSapiens") {
        key_ind <- grepl("All tabula", tabula_key$tabula_sapiens)
    } else{
        key_ind <- grepl(cur_tabula, tabula_key$tabula_sapiens)
    }
    gtex_tissues <- tabula_key$r_names[key_ind]
    genes_of_interest <- unlist(GeneClusterPuller(gtex_tissues, gene_location))
    genes_of_interest <- unique(genes_of_interest)
    ## Filter Seurat data to genes of interest
    genes_of_overlap <- genes_of_interest[genes_of_interest %in% cur_dat@assays$RNA@meta.features$gene_symbol]
    ## Create model matrix
    model_mat <- model.matrix(~ compartment + free_annotation - 1, data = cur_dat@meta.data)
    ## Fix names
    colnames(model_mat) <- gsub(" ", "_", colnames(model_mat))
    colnames(model_mat) <- gsub("free_annotation", "", colnames(model_mat))
    colnames(model_mat) <- gsub("compartment", "", colnames(model_mat))
    ## seperate the matrices on cell types and cell compartments
    compartments <- as.character(unique(cur_dat@meta.data$compartment))
    free_anno <- str_replace_all(as.character(unique(cur_dat@meta.data$free_annotation)),pattern = " ",replacement =  "_")
    model_compartment <- model_mat[,colnames(model_mat) %in% compartments]
    model_celltype <- model_mat[,colnames(model_mat) %in% free_anno]
    ## Run the analysis
    pre_master <- NULL
    ## Don't run all cell types if we're doing the general TabulaSapiens
    if (cur_tabula != "TabulaSapiens") {
        ## For each cell type
        for (col_ind in seq(ncol(model_celltype))) {
            cur_cell <- colnames(model_celltype)[col_ind]
            cell_ind <- model_celltype[,col_ind] == 1
            ## To catch low cell counts
            if (cell_ind < 3) {
                next
            }
            ## Get the names of the cells we're comparing between
            main_cell_names <- rownames(model_celltype)[cell_ind]
            non_cell_names <- rownames(model_celltype)[!cell_ind]
            ## Run ROC analysis
            cell_type_roc <- FindMarkers(cur_dat,
                                      ident.1 = main_cell_names,
                                      only.pos =TRUE,
                                      features = genes_of_overlap,
                                      test.use =  "roc") %>%
                tibble::rownames_to_column("gene_name")
            ## Run wald analysis
            cell_type_wald <- FindMarkers(cur_dat,
                                         ident.1 = main_cell_names,
                                         only.pos =TRUE,
                                         features = genes_of_overlap) %>%
                tibble::rownames_to_column("gene_name")
            ## Combine all the data, adjust and filter
            cur_results <- left_join(cell_type_roc, cell_type_wald) %>%
                mutate(p_val_adj = p.adjust(p_val, method= "bonferroni"),
                       cell_type = cur_cell,
                       tissue = cur_tabula) %>%
                filter(p_val_adj < 0.05)
            ## Add to pre master
            pre_master <- rbind(pre_master, cur_results)
        }
    }
    ## Run comparisons just for cell compartments same code as above
    for (col_ind in seq(ncol(model_compartment))) {
        cur_cell <- colnames(model_compartment)[col_ind]
        cell_ind <- model_compartment[,col_ind] == 1
        if (cell_ind < 3) {
            next
        }
        main_cell_names <- rownames(model_compartment)[cell_ind]
        non_cell_names <- rownames(model_compartment)[!cell_ind]
        cell_type_roc <- FindMarkers(cur_dat,
                                     ident.1 = main_cell_names,
                                     only.pos =TRUE,
                                     features = genes_of_overlap,
                                     test.use =  "roc") %>%
            tibble::rownames_to_column("gene_name")
        cell_type_wald <- FindMarkers(cur_dat,
                                      ident.1 = main_cell_names,
                                      only.pos =TRUE,
                                      features = genes_of_overlap) %>%
            tibble::rownames_to_column("gene_name")
        cur_results <- left_join(cell_type_roc, cell_type_wald) %>%
            mutate(p_val_adj = p.adjust(p_val, method= "bonferroni"),
                   cell_type = cur_cell,
                   tissue = cur_tabula)# %>%
           # filter(p_val_adj < 0.05)
        
        pre_master <- rbind(pre_master, cur_results)
    }
    ## Write out the current Tabula tissue (in case all tabula fails lol)
    write.csv(pre_master,
              file = paste0(output_locations, cur_tabula,"_marker_genes.csv"),
              row.names = FALSE)
    ## Bind to final master table
    master_tables <- rbind(master_tables, pre_master)
}

write.csv(master_tables,
          file = paste0(output_locations, "master_tabula_table_sans_blood_tabula.csv"),
          row.names = FALSE)
