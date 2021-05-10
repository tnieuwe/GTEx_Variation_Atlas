mccall_analysis_corr_step <- function(tissue,
                                      ## Where are the files
                                      path_to_files,
                                      ## What type of variance? Either manual
                                      ## where you select the variance filter
                                      ## using a number, or quantile, which 
                                      ## you give a quantile for a variance 
                                      ## cut off.
                                      variance_type,
                                      ## The value of the cutoff, value's
                                      ## definition changes on the previous
                                      ## argument
                                      variance_value,
                                      correlation = "spearman"){
    ### The purpose of this function is to complete the first and most
    ### computationally intensive section of the McCall pipeline which is
    ### loading in the normalized data, filtering the genes on a variance
    ### threshold and correlating them together.
    
    ## loading the data
    load(paste0(path_to_files,tissue,"-vsd-mean-filtered.rda"))
    vsdMeanFiltered <- generalVSDMeanFiltered
    v <- apply(assay(vsdMeanFiltered),1,var)
    ## Selecting the variance type and creating the index
    if (variance_type == "manual") {
        variance_thresh <- variance_value
        ind <- which(v > variance_thresh)
    } else{
        variance_thresh <- quantile(v, variance_value)
        ind <- which(v > variance_thresh)
        variance_thresh <- round(variance_thresh,2)
    }
    
    ## Creating highvar subset
    vsdHighVar <- assay(vsdMeanFiltered)[ind,]
    rownames(vsdHighVar) <- gtabMeanFiltered$gene_name[ind]
    colnames(vsdHighVar) <- vsdMeanFiltered$SAMPID
    ## Removing duplicated genes
    dim(vsdHighVar)
    vsdHighVar <-  vsdHighVar[!(duplicated(rownames(vsdHighVar))),]
    dim(vsdHighVar)
    vsdHighVarCentered <- vsdHighVar - rowMeans(vsdHighVar)
    
    ## Correlation between high variance genes
    ## This section takes a LONG time if running Kendall's
    gcor <- cor(t(vsdHighVar), method= correlation)
    ## Save data in a list that will be used in all further pipeline steps
    out_list <- list()
    out_list[["correlation_matrix"]] <- gcor
    out_list[["high_var_genes"]] <- vsdHighVar
    out_list[["high_var_centered"]] <- vsdHighVarCentered
    out_list[["tissue"]] <- tissue
    out_list[["correlation"]] <- correlation
    out_list[["var_thresh"]] <- variance_thresh
    return(out_list)
    
}