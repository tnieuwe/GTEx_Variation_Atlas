## Generic functions used throughout the project


#sampids <- stab$SAMPID

GTEx_SAMPID_to_SUBJID <- function(sampids){
    ## A simple function to quickly turn SAMPIDs into SUBJIDs to connect
    ## individuals to their phenotype
  stringr::str_extract(string = sampids, pattern = "GTEX-[[:alnum:]]*")
}


ZscoreProfilePuller <- function(r_name_vect, location){
    ### This function is used to pull Z-score of different tissues easily
    profile_list <- list()
    for (tiss in r_name_vect) {
        cur_frame <- read.csv(paste0(location,
                                     "kendall-",
                                     tiss,
                                     "-cluster-profiles.csv")) %>%
            mutate(SUBJID = GTEx_SAMPID_to_SUBJID(X)) %>%
            select(SUBJID, SAMPID = "X", everything())
        profile_list[[tiss]] <- cur_frame
    }
    return(profile_list)
}

GeneClusterPuller <- function(r_name_vect, location){
    ### This function is used to pull the genes of tissue clusters
    tiss_list <- list()
    for (tiss in r_name_vect) {
        cur_frame <- read.csv(paste0(location,
                                     "kendall-",
                                     tiss,
                                     "-gene-clusters-high-variance.csv"))
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

AllGeneClusterPuller <- function(location){
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

ClusterGeneOverlaps <- function(genes_of_int,
                                all_clusters,
                                only_hits = TRUE){
  ### This function finds what clusters overlap with a provided list of genes
  res_mat <- matrix(nrow = sum(unlist(lapply(all_clusters, length))), ncol = 2)
  rownames(res_mat) <- seq(nrow(res_mat))
  overall_ind = 1
  for (tissue_ind in seq_along(all_clusters)) {
    tiss_name = names(all_clusters)[tissue_ind]
    cur_tiss <- all_clusters[[tissue_ind]]
    for (cluster_ind in seq_along(cur_tiss)) {
      clust_name = names(cur_tiss)[cluster_ind]
      cur_cluster <- cur_tiss[[cluster_ind]]
      overlapping_genes <- intersect(cur_cluster, genes_of_int)
      n_of_overlap <- length(overlapping_genes)
      ## Naming
      rownames(res_mat)[overall_ind] = paste(tiss_name, clust_name,sep = "_")
      res_mat[overall_ind, 1] = n_of_overlap
      if (length(overlapping_genes) == 0) {
        overlapping_genes = ""
      }
      res_mat[overall_ind, 2] = paste(overlapping_genes, collapse = "; ")
      
      overall_ind = overall_ind + 1
    }
  }
  ## Clean up table
  res_df <- as.data.frame(res_mat)
  res_df$V1 <- as.numeric(res_df$V1)
  colnames(res_df) <- c("n_intersect", "found_genes")
  
  if (only_hits == TRUE) {
    res_df <- dplyr::filter(res_df, n_intersect > 0)    
  }
  return(res_df)
}

CommonGenes <- function(gene_of_int, cluster_list, x_clust = TRUE){
  ### The purpose of this function is to find out what genes frequently
  ### show up with another gene of interest. We also have the functionality
  ### to ignore X clusters as those aren't "real clusters". 
  final_vect <- c()
  ## using a for loop 
  for (tiss_ind in seq_along(cluster_list)) {
    cur_tiss <- cluster_list[[tiss_ind]]
    keep_ind <- lapply(cur_tiss, function(x){gene_of_int %in% x})
    if (x_clust == FALSE & keep_ind[['X']] == TRUE ){next}
    
    final_vect<- append(final_vect, unlist(cur_tiss[unlist(keep_ind)]))
  }
  pre_final_table <- as.data.frame(table(final_vect))
  pre_final_table <- pre_final_table[pre_final_table$Freq!=1,]
  final_table<- pre_final_table[order(pre_final_table$Freq,decreasing = T),]
  return(final_table)
}

RemoveSemicolon <- function(gene_list){
  ### The purpose of this code is to make it easier to copy and paste
  ### gene lists into tools such as HPAStainR and GO Ontology
  cat(paste0(unlist(strsplit(gene_list,split =  ";")), "\n"))
}


GenesFromTissues <- function(tiss_of_int, genes_of_int, key_dat, vst_location,
                             run_all_available = FALSE){
  ### The purpose of this function is to easily load in normalized RNA-seq
  ### data to compare gene expression. The main idea is to see if the
  ### expression correlates across samples.
  ## load required packages
  require(DESeq2)
  require(tidyverse)
  ## Filter tissue based on any of the columns, this allows us to use abbreviations
  ## r names or names proper
  if (run_all_available == TRUE){
    tiss_of_int <- c()
    for (tissue in unique(key_dat$r_names)) {
      tiss_path <- paste0(vst_location,tissue,"-vsd-mean-filtered.rda")
      if (file.exists(tiss_path)) {
        tiss_of_int <- append(tiss_of_int, tissue)
      }
    }
    
  }
  cur_keys <- key_dat %>% select(abrreviated, r_names, official_name) %>%
    unique() %>%
    filter(case_when(
      abrreviated %in% tiss_of_int |
        r_names %in% tiss_of_int |
        official_name %in% tiss_of_int ~ TRUE,
      TRUE ~ FALSE
    ))
  cur_keys <- cur_keys[order(tiss_of_int, cur_keys$r_names),]
  
  ## Prepare returned data frame
  final_gene_dat <- NULL
  for (tiss_ind in seq_along(tiss_of_int)) {
    ## Select and load tissues which exist
    tissue <- cur_keys$r_names[tiss_ind]
    abbrev <- cur_keys$abrreviated[tiss_ind]
    tiss_path <- paste0(vst_location,tissue,"-vsd-mean-filtered.rda")
    ## Catch missing vsts
    if (!file.exists(tiss_path)) {
      stop(paste0(tissue, " PATH does not exist"))
    }
    ## Load data and do a quick sanity check
    load(tiss_path)
    if(!all(gtabMeanFiltered$gene_id == rownames(assay(generalVSDMeanFiltered)))){
      stop(paste0(tissue, " fails sanity check"))
    }
    ## Index on genes of interest
    ind <- gtabMeanFiltered$gene_name %in% genes_of_int
    final_gtab <- gtabMeanFiltered[ind,]
    ## Take out gene data and name it
    gene_dat <- t(assay(generalVSDMeanFiltered)[ind,,drop =FALSE])
    colnames(gene_dat) <- final_gtab$gene_name
    ## If more than one tissue bind the tissue name
    if (length(tiss_of_int) > 1) {
      gene_dat <- cbind(as.data.frame(gene_dat), tissue)
      ## if first tissue just place it in final object
      if (is.null(final_gene_dat)) {
        final_gene_dat <- gene_dat
        next
      }
      #Use bind rows from dplyr in case not all genes are shareed
      final_gene_dat <- bind_rows(final_gene_dat, gene_dat) %>%
        select(everything(), tissue)
    } else{
      colnames(gene_dat) <-final_gtab$gene_name
      gene_dat <- rownames_to_column(as.data.frame(gene_dat), var = "SAMPID")
      final_gene_dat <- gene_dat
    }
    
  }
  if (!("SAMPID" %in% colnames(final_gene_dat))) {
    final_gene_dat <- final_gene_dat %>% rownames_to_column(var = "SAMPID")
  }
  ## 
  final_gene_dat <- final_gene_dat  %>%
    mutate(SUBJID = GTEx_SAMPID_to_SUBJID(SAMPID)) %>%
    select(SUBJID, SAMPID, everything())

  
  
  return(final_gene_dat)
  
}


GenesFromTissWiden <- function(GenesFromTissues_output){
  pre_bound <- GenesFromTissues_output
  all_bound <- c()
  for (cur_ind in seq_along(unique(pre_bound$tissue))) {
    cur_tiss <- unique(pre_bound$tissue)[cur_ind]
    temp_dat <- filter(pre_bound, tissue == cur_tiss) %>%
      select(-tissue)
    new_subj <- GTEx_SAMPID_to_SUBJID(temp_dat$SAMPID)
    colnames(temp_dat) <- paste0(colnames(temp_dat), "_",cur_tiss)
    temp_dat <- cbind(new_subj, temp_dat)
    if (cur_ind == 1) {
      all_bound <- temp_dat
      next
    } else{
      all_bound <- left_join(all_bound, temp_dat)
    }
    
  }
  return(all_bound)
}
