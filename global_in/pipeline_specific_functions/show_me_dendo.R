show_me_dendo <-function(cluster_input,
                         ## What cluster do you want to see?
                         cluster_select,
                         ## Where do you want the line?
                         draw_line = NULL){
    ### The purpose of this function is to see dendrogram of clusters so you
    ### can manually break them up with the function `clust_breaker()`. This
    ### function allows you to draw a line on the dendrogram to see where it
    ### would cut. Uses list input generated from McCall pipeline functions. 
    gene_lists <- cluster_input$cluster_genes
    gene_clust <- gene_lists[,stringr::str_detect(colnames(gene_lists),
                                                  cluster_select)]
    gene_clust <- as.vector(gene_clust)
    cleaned_genes <- gene_clust[gene_clust != ""]
    cor_mat <- cluster_input$correlation_matrix
    sub_cor <- cor_mat[cleaned_genes, cleaned_genes]
    dendo_sub <- hclust(d=as.dist(1-abs(sub_cor)), method = "average")
    if (is.null(draw_line)) {
        plot(dendo_sub, h = -1,cex = 0.6)
        
    } else{
        plot(dendo_sub, h = -1,cex = 0.6)
        abline(h = draw_line)
    }
    
}