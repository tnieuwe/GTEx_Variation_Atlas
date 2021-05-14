clust_breaker <- function(cluster_input,
                          cluster_select,
                          ## Select if you wnat to split on N of clusters or
                          ## height.
                          number_o_height,
                          ## The N of cluster or height of cut you want
                          value,
                          ## How many genes you want in each cluster
                          min_clust = 6){
    ### The purpose of this code is to break up clusters manually that weren't
    ### broken by previous methods. These breaks are completely arbitrary, but
    ### we can find interesting variation inside of them. Eg. Break cluster B
    ### from the kidney results in proximal convoluted tubule and podocyte
    ### clusters. Requires list inputs created by McCall pipeline.
    
    ## Load in data and cut dendo
    gene_lists <- cluster_input$cluster_genes
    gene_clust <- gene_lists[,stringr::str_detect(colnames(gene_lists), cluster_select)]
    gene_clust <- as.vector(gene_clust)
    cleaned_genes <- gene_clust[gene_clust != ""]
    cor_mat <- cluster_input$correlation_matrix
    ## Make sub correlation matrix
    sub_cor <- cor_mat[cleaned_genes, cleaned_genes]
    dendo_sub <- hclust(d=as.dist(1-abs(sub_cor)), method = "average")
    ## Here we cut the tree as we prefer
    if (number_o_height == "height") {
        new_clust <- cutree(dendo_sub, h = value)
        
    }else{
        new_clust <- cutree(dendo_sub, k = value)
    }
    ## Filter clusters and make gclust table
    clust_to_keep <- as.numeric(names(table(new_clust)[table(new_clust) >= min_clust]))
    clust_filt <- new_clust[new_clust %in% clust_to_keep]
    clust_table <- as.data.frame(clust_filt)
    gclusts <- matrix("", nrow=max(table(clust_filt)), ncol=length(table(clust_filt)))
    for(k in 1:length(table(clust_filt))){
        ## The index below is broken I will fix it 5/4/2021
        cur_genes <- rownames(clust_table)[clust_table == k]
        #gclusts[1:sum(cnames==tmp), k] <- rownames(gcor)[cnames==tmp]
        gclusts[,k] <- c(cur_genes, rep("", nrow(gclusts) - length(cur_genes)))
        
    }
    ## Make Z-scores for new clusts
    vsdHighVarCenteredClust <- cluster_input$high_var_clust_cent
    ## extract expression for each cluster, convert to zscores, and summarize profile
    clusterProfiles <- matrix(nrow=ncol(vsdHighVarCenteredClust), ncol=(ncol(gclusts)))
    for(k in 1:ncol(clusterProfiles)){
        # ind <- which(cnamesClust == LETTERS[k])
        
        #New code picks from gclusts
        ind <- rownames(vsdHighVarCenteredClust) %in%  gclusts[,k]
        
        #For normal clusters 
        if (length(unique(gclusts[,k])) != 2) {
            
            etmp <- vsdHighVarCenteredClust[ind,]
            ztmp <- etmp / apply(etmp, 1, mad)
            #Below removes inf values
            ztmp <- ztmp[is.finite(rowSums(ztmp)),]
            
            clusterProfiles[,k] <- colMeans(ztmp)
        } else {
            #For single genes such as XIST
            etmp <- vsdHighVarCenteredClust[ind,]
            clusterProfiles[,k] <- etmp
            
        }
        
    }
    
    new_col_names <- paste0(rep(cluster_select, ncol(gclusts)), seq(ncol(gclusts)))
    rownames(clusterProfiles) <- colnames(vsdHighVarCenteredClust)
    colnames(gclusts) <- new_col_names
    colnames(clusterProfiles) <- new_col_names
    ## Return different list output when compared to input list
    dat_out <- list()
    dat_out[["sub_genes"]] <- gclusts
    dat_out[["sub_profile"]] <- clusterProfiles
    dat_out[["dendo_sub"]] <- dendo_sub
    dat_out[["sub_cor"]] <- sub_cor
    return(dat_out)
}