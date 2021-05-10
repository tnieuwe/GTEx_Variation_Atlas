mccall_profile_maker <- function(sub_clust_out){
    ### This code generates for each cluster a Z-score using the normalized
    ### gene data. These Z-scores are called 'cluster profiles' in much of this
    ### analysis and likely paper. 
    gclusts <- sub_clust_out$cluster_genes
    vsdHighVarCenteredClust <- sub_clust_out$high_var_clust_cent
    ## extract expression for each cluster, convert to zscores, and summarize profile
    clusterProfiles <- matrix(nrow=ncol(vsdHighVarCenteredClust), ncol=(ncol(gclusts) - 1))
    for(k in 1:ncol(clusterProfiles)){
        # ind <- which(cnamesClust == LETTERS[k])
        
        #New code picks from gclusts
        ind <- rownames(vsdHighVarCenteredClust) %in%  gclusts[,1 + k]
        
        #FOr normal clusters 
        if (length(unique(gclusts[,1 + k])) != 2) {
            
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
    rownames(clusterProfiles) <- colnames(vsdHighVarCenteredClust)
    colnames(clusterProfiles) <- colnames((gclusts))[-1]
    ## Save and return output
    sub_clust_out[["sample_profiles"]] <- clusterProfiles
    sub_clust_out
}