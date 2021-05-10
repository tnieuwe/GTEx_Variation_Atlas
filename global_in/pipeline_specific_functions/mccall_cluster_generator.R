mccall_cluster_generator <- function(cor_out,
                                     ## Cluster breaker if set to TRUE includes
                                     ## the agglomerative approach to break up
                                     ## subclusters.
                                     cluster_breaker = FALSE,
                                     required_genes = 6){
    ### This function loads in the output of mccall_analysis_corr_step and
    ### creates clusters based on a tau critical value generated based off the
    ### number of samples and the number of genes in the dataset. Also includes
    ### an agglomerative approach to break up subclusters based on a quantile
    ### correlation threshold specifc to each subcluster.
    
    ## Load in data from input list
    vsdHighVar <- cor_out[["high_var_genes"]]
    gcor <- cor_out[['correlation_matrix']]
    ## Creat tau crit
    ngene <- nrow(vsdHighVar)
    nsamp <- ncol(vsdHighVar)
    tauCrit <- round(qKendall(0.01/((ngene^2)-ngene)/2, nsamp, lower.tail = F), 2)
    
    ## cluster based on absolute correlation distance
    dendo <- hclust(d=as.dist(1-abs(gcor)), method = "average")
    
    ## clean up heatmap to focus on interesting gene clusters
    clusts <- cutree(dendo, h=1-tauCrit)
    tabclusts <- table(clusts)
    ind <- which(tabclusts >= required_genes) #arbitrary
    clusts[!(clusts %in% names(tabclusts)[ind])] <- 0
    ntmp <- names(clusts)
    clusts <- as.numeric(as.factor(clusts))-1
    names(clusts) <- ntmp
    
    ### BEGINNING OF CLUSTER BREAKER ###
    if (cluster_breaker == TRUE) {
        for (each_clust in sort(unique(clusts))[-1]) {
            ## Subset out clusters to test each of them
            genes <- names(clusts[clusts %in% each_clust])
            sub_clust <- (gcor[genes,genes])
            cur_clust_mean <- mean(abs(sub_clust))
            ## Create a threshold that determine if we slice a cluster or not
            break_down_threshold <- quantile(abs(gcor), .85)
            if (cur_clust_mean < break_down_threshold) {
                ## Create cluster  dendrogram used to split cluster into
                ## subclusters
                sub_dendo <- hclust(d=as.dist(1-abs(sub_clust)), method = "average")
                cur_thresh <- quantile(abs(sub_clust), .85) ## Quantile value is arb
                cur_min_corr <- 1
                starting_cut <- 0
                ## This loops starts at the height of 0 on the dendrogram where
                ## every cluster is one gene with a correlation of 1. It then 
                ## moves up in steps of 0.001 and recalculates the correlation
                ## of all generated groups. As long as the minimum correlations
                ## of the groups remain below the threshold created above, the
                ## line will continue to move. Once the threshold is passed
                ## the loops stop and those are our final clusters.
                while (cur_thresh <= cur_min_corr) {
                    clusts_sub <- cutree(sub_dendo, h = starting_cut)
                    corr_vals <- sapply(unique(clusts_sub), FUN = function(x){
                        cur_clust <- clusts_sub[clusts_sub == x]
                        cur_genes <- names(cur_clust)
                        mean(abs(sub_clust[cur_genes,cur_genes]))
                    })
                    cur_min_corr = min(corr_vals)
                    ## Arb step, could do some gradient descent stuff here if 
                    ## runtime ever becomes prohibitive
                    starting_cut = starting_cut + 0.001
                }
                all_clusts <- table(clusts_sub)
                ## arbitrary N of required genes for a clusters
                ind <- all_clusts >= required_genes
                clusts_filt <- all_clusts[ind]
                clusts_sub[!(clusts_sub %in% names(clusts_filt))] <- 0
                ## Remove current clust from clusts
                clusts <- clusts[!(clusts == each_clust)]
                ## Split non-clust and clustered genes
                no_clust_only <- clusts_sub[clusts_sub == 0]
                true_clust_only <- clusts_sub[clusts_sub != 0]
                ## Generate new cluster numbers based on the the current cluster
                ## numbers.
                new_nums <- seq(max(clusts) + 1,
                                max(clusts) + length(unique(true_clust_only)))
                ## Create dummy variable to prevent folding in of clusters
                dummy_true <- paste0(true_clust_only, "A")
                names(dummy_true) <-names(true_clust_only)
                for (ind in seq(length(unique(true_clust_only)))) {
                    dummy_true[true_clust_only == sort(unique(true_clust_only))[ind]] <- new_nums[ind]
                }
                numeric_true <- as.numeric(dummy_true)
                names(numeric_true) <- names(dummy_true)
                clusts_sub <- c(no_clust_only, numeric_true)
                ## Rejoin non-clustering genes
                clusts <- c(clusts, clusts_sub)
            }
            
        }
        ## Rename clusters here
        new_nums <- seq(0, length(unique(clusts)) - 1)
        for (ind in seq(length(unique(clusts)))) {
            clusts[clusts == unique(clusts)[ind]] <- new_nums[ind]
        }
    }
    ### END OF CLUSTER BREAKER ###
    #Naming clusters based on their clust values
    n <- length(table(clusts))
    cnames <- rep("", length(clusts))
    ## Replace this with an array based approach
    for (clst in 1:(n)){
        cnames[clusts==clst] <- LETTERS[clst]
    }
    
    anno = data.frame("cluster"=factor(cnames))
    rownames(anno) = names(clusts)
    ## Save data to same list and return said list
    out_list <- cor_out
    out_list[["anno_df"]] <- anno
    out_list[["cluster_names"]] <- cnames
    out_list[["clust_num"]] <- n
    out_list[["dendogram"]] <- dendo
    
    return(out_list)
}