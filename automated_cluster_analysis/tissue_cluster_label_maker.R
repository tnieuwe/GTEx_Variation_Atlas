### The purpose of this file is to create labels for each tissue and cluster
### for easier data querying in the future.

## Load packages
library(stringr)
library(dplyr)
library(tibble)
library(tidyr)
## load data

## Tissue dictionary
tiss_dict <- read.csv("automated_cluster_analysis/data_in/ordered_all_tiss_mat_names.csv")

## New clusters
tiss_list <- read.csv(file = "global_in/general_list_test_v8.txt",header = F)[,1]
drop_tiss <- c("cervix_ectocervix", "cervix_endocervix", "bladder","fallopian_tube")
tiss_list <- tiss_list[!(tiss_list %in% drop_tiss)]

list_dat <- sapply(tiss_list, function(x){
    ## Go through list to load in files
    load_file <-  paste0("MARCC output/variance_genes_kendall/kendall-",
                         x,
                         "-gene-clusters-high-variance.csv")
    new_dat <- read.csv(load_file)
    all_col_names <- colnames(new_dat)
    # for (cur_col_name in all_col_names) {
    #     ## Below code collapses the negatively correlating gene clusters together
    #     #cur_col_name <- names(new_dat)[col_ind]
    #     # if (is.na(cur_col_name)) {
    #     #     break
    #     # }
    #     if (str_detect(cur_col_name, "2")) {
    #         next
    #     }
    #     if(str_detect(cur_col_name, "1")){
    #         flatten_col_name <- str_replace(cur_col_name, ".1", "")
    #         col_1 <- new_dat[[paste0(flatten_col_name, ".1")]]
    #         col_2 <- new_dat[[paste0(flatten_col_name, ".2")]]
    #         join_col <- append(col_1, col_2)
    #         join_col <- join_col[join_col != ""]
    #         join_col <- append(join_col, rep("", length(col_1) - length(join_col)))
    #         remove_old_ind <- !startsWith(colnames(new_dat), flatten_col_name)
    #         new_dat <- new_dat[remove_old_ind]
    #         new_dat <- cbind(new_dat, join_col)
    #         names(new_dat)[ncol(new_dat)] <- flatten_col_name
    #         
    #     }
    # }
    ## Reorder colnames
    pre_order <- order(colnames(new_dat))
    order_ind <- pre_order[-length(pre_order)]
    order_ind <- append(1, order_ind)
    new_dat <- new_dat[,order_ind]
    df_list <- NULL
    unclustered_vect <- c()
    
    for (col_ind in seq(ncol(new_dat))) {
        new_col <- new_dat[,col_ind]
        ## Remove non-clustering clusters
        # if (colnames(new_dat)[col_ind] != "X") {
        #   df_list[[colnames(new_dat)[col_ind]]] <- new_col[new_col != ""]  
        # }
        ## Keep non-clustering genes and
        df_list[[colnames(new_dat)[col_ind]]] <- new_col[new_col != ""]  
    }
    return_list <- NULL
    return_list[[x]] <- df_list
    
})
### Make labels
tiss_dict <- tiss_dict[tiss_dict$r_names %in% names(list_dat),]

pre_df <- lapply(list_dat, names)
key_table <- NULL
for (ind in seq_along(pre_df)) {
    true_ind = 1
    cur_tiss <- names(pre_df)[ind]
    temp_table <- NULL
    for (jnd in seq(length(pre_df[[ind]]))) {
        cur_clust = pre_df[[ind]][[jnd]]
        new_row <- cbind(cur_tiss, cur_clust)
        key_table<- rbind(key_table, new_row)
        
    }
    true_ind = true_ind + 1
}

key_table <- key_table %>% as.data.frame() %>%
    rename(r_names = "cur_tiss", cluster = "cur_clust") %>%
    left_join(tiss_dict)
final_key_table <- key_table %>%
    mutate(label = str_replace(paste(abrreviated, cluster,sep  = "_"), "[.]", "_")) %>%
    select(label, abrreviated, cluster, r_names, official_name)
## Write the key table output
write.csv(final_key_table,
          file = "automated_cluster_analysis/data_out/gtex_tiss_clust_key_table.csv",row.names = F)

## Make a final output including the genes pasted together
genes_per_group <- lapply(list_dat, function(x){
    lapply(x, function(y){
        paste(y,collapse = ";")
    })
})
gene_group_df <- (as.data.frame((unlist(genes_per_group)))) %>%
    rownames_to_column("label") %>%
    rename(genes = "(unlist(genes_per_group))") %>%
    mutate(label = str_replace(label, "[.]", ";")) %>%
    separate(col = "label", into = c("r_names", "cluster"), sep = ";") %>%
    left_join(final_key_table, ., by = c("r_names", "cluster")) %>%
    select("label", "genes")
write.csv(gene_group_df,
            file = "automated_cluster_analysis/data_out/labels_and_genes.csv",
            row.names = FALSE,quote = FALSE)    
