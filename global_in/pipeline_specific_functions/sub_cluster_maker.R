sub_cluster_maker <- function(clust_initial){
    ### This function splits the clusters created by absolute values into two
    ### seperate cluster eg. X-1 and X-2. This is to make the clusters more
    ### informative. 
    ## Load in the data
    cnames <- clust_initial$cluster_names
    anno <- clust_initial$anno_df
    gcor <- clust_initial$correlation_matrix
    vsdHighVar <- clust_initial[["high_var_genes"]]
    high_var_centered <- clust_initial$high_var_centered
    ## Generate gclusts matrix, which includes the genes and their respective
    ## clusters.
    gclusts <- matrix("", nrow=max(table(cnames)), ncol=length(table(cnames)))
    colnames(gclusts) <- levels(anno$cluster)
    for(k in 1:length(table(cnames))){
        tmp <- colnames(gclusts)[k]
        cur_genes <- rownames(anno)[anno$cluster == tmp]
        #gclusts[1:sum(cnames==tmp), k] <- rownames(gcor)[cnames==tmp]
        gclusts[1:sum(cnames==tmp), k] <- cur_genes
        
    }
    
    ## Below code creates subclusters, splitting on the absolute values
    for (column in LETTERS[1:(length(unique(cnames)) -1 )]) {
        #Pull out the first column
        test_col <- gclusts[,column]
        first_gene <- test_col[1]
        #Make an index to subset gcor to just a single clusters
        ind <- which(rownames(gcor) %in% test_col)
        sub_cor <- gcor[ind,ind]
        #Test if the cluster has any negatively correlating genes
        ## This doesn't work let me try a new method generating a new dendo
        ind2 <- sub_cor[,first_gene] < 0
        #If there are negatively correlating genes...
        if (sum(ind2) > 0) {
            dendo_sub <- hclust(d=as.dist(1-(sub_cor)), method = "average")
            
            ## clean up heatmap to focus on interesting gene clusters
            sub_clusts <- cutree(dendo_sub,k = 2)
            
            #generate sub cluster 1
            sub_1 <- names(sub_clusts)[sub_clusts == 1]
            #Subcluster 2
            sub_2 <- names(sub_clusts)[sub_clusts == 2]
            ### Triple Check by seeing if output pos or neg correlates
            ### Added 12/21/2021 late in the process to catch strange cluster
            ### splitting
            temp_z <- cbind(colMeans(high_var_centered[sub_1,,drop = F]),
                            colMeans(high_var_centered[sub_2,,drop = F]))
            do_they_correlate <- cor(temp_z)[1,2] > 0 
            
            ###
            if (do_they_correlate == FALSE) {
                #Remover old column to replace with new columns
                gclusts <- as.matrix(gclusts[,!(colnames(gclusts) == column)])
                #Make a single matrix of same rows as gclusts, one for each sub-cluster. 
                matrix_1  <- as.matrix(c(sub_1 ,rep("", nrow(gclusts) - length(sub_1))), )
                colnames(matrix_1) <- paste0(column,"-1")
                matrix_2  <- as.matrix(c(sub_2 ,rep("", nrow(gclusts) - length(sub_2))), )
                colnames(matrix_2) <- paste0(column,"-2")
                #Add subclusters to gclusts
                gclusts <-   cbind(gclusts, matrix_1, matrix_2)
            }
            
        }
        
    }
    
    #Realphebatize columns for sake of beauty 
    gclusts <- gclusts[,order(colnames(gclusts))]
    
    ## Stratify cluster inclusion
    ind <- which(rownames(vsdHighVar) %in% rownames(anno)[anno$cluster != ""])
    vsdHighVarClust <- clust_initial$high_var_genes[ind,]
    vsdHighVarCenteredClust <- clust_initial$high_var_centered[ind,]
    vsdHighVarNoise <- clust_initial$high_var_genes[-ind,]
    
    ## Return output to the input list
    clust_initial[["cluster_genes"]] <- gclusts
    clust_initial[["high_var_clusters"]] <- vsdHighVarClust
    clust_initial[["high_var_clust_cent"]] <- vsdHighVarCenteredClust
    clust_initial[["high_var_noise"]] <- vsdHighVarNoise
    clust_initial
}