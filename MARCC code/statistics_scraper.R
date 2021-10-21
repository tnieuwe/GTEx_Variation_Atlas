##### Statistic Scraper #####
### The purpose of this script is to scrape the required statistics
### from MARCC for the paper this includes

### Load functions and libararies
library(DESeq2)

AllGeneClusterPuller <- function(location, x_clust = TRUE){
    ### This function is used to pull the genes of ALL tissue clusters
    r_name_vect <- c("adipose_subcutaneous", "adipose_visceral__omentum", "adrenal_gland", 
                     "artery_aorta", "artery_coronary", "artery_tibial", "brain_amygdala", 
                     "brain_anterior_cingulate_cortex__ba24", "brain_caudate__basal_ganglia", 
                     "brain_cerebellar_hemisphere", "brain_cerebellum", "brain_cortex", 
                     "brain_frontal_cortex__ba9", "brain_hippocampus", "brain_hypothalamus", 
                     "brain_nucleus_accumbens__basal_ganglia", "brain_putamen__basal_ganglia", 
                     "brain_spinal_cord__cervical_c_1", "brain_substantia_nigra", 
                     "breast_mammary_tissue", "cells_ebv_transformed_lymphocytes", 
                     "cells_cultured_fibroblasts", "colon_sigmoid", "colon_transverse", 
                     "esophagus_gastroesophageal_junction", "esophagus_mucosa", "esophagus_muscularis", 
                     "heart_atrial_appendage", "heart_left_ventricle", "kidney_cortex", 
                     "liver", "lung", "minor_salivary_gland", "muscle_skeletal", "nerve_tibial", 
                     "ovary", "pancreas", "pituitary", "prostate", "skin_not_sun_exposed__suprapubic", 
                     "skin_sun_exposed__lower_leg", "small_intestine_terminal_ileum", 
                     "spleen", "stomach", "testis", "thyroid", "uterus", "vagina", 
                     "whole_blood")
    tiss_list <- list()
    for (tiss in r_name_vect) {
        cur_frame <- read.csv(paste0(location,
                                     "kendall-",
                                     tiss,
                                     "-gene-clusters-high-variance.csv"))
        ## Add X cluster removal option
        if (x_clust == FALSE) {
            x_ind <-colnames(cur_frame) %in% "X" 
            cur_frame <- cur_frame[,!x_ind]
        }
        
        clust_list <- list()
        for (cur_col in seq(ncol(cur_frame))) {
            gene_column <- cur_frame[,cur_col]
            gene_column <- gene_column[gene_column !=""]
            column_name <- colnames(cur_frame)[cur_col]
            clust_list[[column_name]] <- gene_column
        }
        tiss_list[[tiss]] <- clust_list
    }
    return(tiss_list)
}



### N of genes at each stage of analysis
cluster_locations <- "~/work2/tnieuwe1/data/gtex_v8/kendall_gene_clusters/"
cluster_list <- AllGeneClusterPuller(location = cluster_locations)
cluster_list_no_x <- AllGeneClusterPuller(location = cluster_locations,x_clust = F)
all_tissues <- read.table("~/work2/tnieuwe1/data/gtex_v8/gen_tiss_lists/general_list_test_v8.txt")[,1]
data_out <- NULL
for (tiss_ind in seq_along(all_tissues)) {
    cur_tiss <- all_tissues[tiss_ind]
    mean_filtered_location <- paste0("~/work2/tnieuwe1/data/gtex_v8/diff_samps/",
                                     cur_tiss, "/", cur_tiss,
                                     "-vsd-mean-filtered.rda")
    load(mean_filtered_location)
    ## Move the object into a "general" file and get variance
    vsdMeanFiltered <- get(paste0(cur_tiss,"VSDMeanFiltered"))
    v <- apply(assay(vsdMeanFiltered),1,var)
    variance_thresh <- quantile(v, .98)
    ## N of filtered genes
    mean_filtered_genes <- nrow(gtabMeanFiltered)
    ## Var filtered genes
    all_var_genes <- length(unlist(cluster_list[[cur_tiss]]))
    ## Cluster genes lsit
    all_clustered_genes <- length(unlist(cluster_list_no_x[[cur_tiss]]))
    new_row <- cbind(cur_tiss, mean_filtered_genes,
                     all_var_genes, all_clustered_genes, variance_thresh)
    data_out <- rbind(data_out, new_row)
}
colnames(data_out) <- c("tissue", "mean_filtered_genes",
                        "variance_filtered_genes", "all_clustered_genes",
                        "variance_threshold")

write.csv(data_out, "~/work2/tnieuwe1/data/gtex_v8/statistics_of_clusters.csv",row.names = F)