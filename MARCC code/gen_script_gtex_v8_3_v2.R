### general analysis

## Load libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(SuppDists)
library(reshape2)
## Below are added for umap
library(umap)
library(dplyr)
library(ggplot2)
library(tibble)
library(stringr)
## Load function
source("~/work2/tnieuwe1/data/gtex_v8/pipeline_scripts/mccall_analysis_corr_step.R")
source("~/work2/tnieuwe1/data/gtex_v8/pipeline_scripts/mccall_cluster_generator.R")
source("~/work2/tnieuwe1/data/gtex_v8/pipeline_scripts/mccall_profile_maker.R")
source("~/work2/tnieuwe1/data/gtex_v8/pipeline_scripts/sub_cluster_maker.R")
source("~/work2/tnieuwe1/data/gtex_v8/pipeline_scripts/sub_cluster_maker.R")
source("~/work2/tnieuwe1/data/gtex_v8/pipeline_scripts/unfiltered_heatmap.R")

## Set up strings
correlation_type = "kendall"
profiles_path <- paste0("~/work2/tnieuwe1/data/gtex_v8/diff_samps/general/",
                        correlation_type, "-general-cluster-profiles.csv")
gene_clust_path <- paste0("~/work2/tnieuwe1/data/gtex_v8/diff_samps/general/",
                          correlation_type, "-general-gene-clusters-high-variance.csv")
heatmap_path <- paste0("~/work2/tnieuwe1/data/gtex_v8/diff_samps/general/",
                       correlation_type, "-general-between-gene-correlation","
                       -high-variance-genes.pdf")
pipeline_rda_path <-  paste0("~/work2/tnieuwe1/data/gtex_v8/pipeline_rdas/",
                             correlation_type, "-general-pipeline-result",".rda")
umap_path <-  paste0("~/work2/tnieuwe1/data/gtex_v8/umaps_out/",
                             correlation_type, "-general-umap-result",".pdf")


## Load data
# load(file="gene_analysis_comparisons/data_in/general-vsd-mean-filtered.rda")
# vsdMeanFiltered <- generalVSDMeanFiltered
pipeline_dat <- mccall_analysis_corr_step("general",
                                     path_to_files = "~/work2/tnieuwe1/data/gtex_v8/diff_samps/general/",
                                     variance_type = "quantile",
                                     variance_value = .98,
                                     correlation = correlation_type,
                                     MARCC = TRUE)
pipeline_dat <- mccall_cluster_generator(pipeline_dat,
                                         cluster_breaker = TRUE,
                                         sub_clust_cut_off = .8,
                                         agglo_cut_off = .75)
pipeline_dat <- sub_cluster_maker(pipeline_dat)
pipeline_dat <- mccall_profile_maker(pipeline_dat)

## Make figures and save data
unfiltered_heatmap(pipeline_dat,
                   save_file = heatmap_path,
                   row_col_font = 2) 
write.csv(pipeline_dat[["sample_profiles"]], file=profiles_path, quote=FALSE)
write.csv(pipeline_dat[["cluster_genes"]], file = gene_clust_path,
          row.names = FALSE)
save(pipeline_dat, file = pipeline_rda_path)

### Add umap creation ###
## Set seed for umap 
set.seed(123)

## Load in gnenes we're dropping
umap_drop <- read.table(file = "~/work2/tnieuwe1/data/gtex_v8/umap_remove_genes.tsv",
           header = TRUE, sep = "\t")
if ("general" == "pancrease") {
    umap_drop<- umap_drop[,-2,drop=FALSE]
}
pre_filt_drop <- unlist(umap_drop)
filt_drop <- pre_filt_drop[pre_filt_drop != ""]
## Get columns to remove and remove them
ind_remove <- !as.logical(colSums(apply(pipeline_dat$cluster_genes,MARGIN = 2,
                                        FUN = function(x){x %in% filt_drop})))
unclean_gene_list <- as.vector(pipeline_dat$cluster_genes[,ind_remove])
clean_gene_list <- unclean_gene_list[unclean_gene_list != ""]
## Filter vst data to only these clusters
unknown_clust_genes <- pipeline_dat$high_var_centered[rownames(pipeline_dat$high_var_centered) %in% clean_gene_list,]
## Make maleable profile df
profile_df <- as.data.frame(pipeline_dat$sample_profiles) %>%
    rownames_to_column("sample")
## Make upmap
norm_umap <- umap(t(unknown_clust_genes))
## Pull matrix
umap_matrix <- norm_umap$layout
## Prepare for ggplot2
umap_for_plot   <- as.data.frame(umap_matrix) %>% rownames_to_column(var = "sample") %>%
    left_join(profile_df)
## Some renaming for easier looping
colnames(umap_for_plot)[2:3] <- c("umap_1", "umap_2")
colnames(umap_for_plot) <- str_replace(colnames(umap_for_plot), pattern = "-", ".")
clusters <- colnames(umap_for_plot)[4:ncol(umap_for_plot)]
## Loop data
plot_list <- list()
for (clusts in clusters) {
    plt <-  ggplot(umap_for_plot, aes_string(x= "umap_1", y = "umap_2", color = clusts)) +
        geom_point() +
        theme_classic()
    plot_list[[clusts]] <- plt
}
## Save data
pdf(umap_path, width = 6, height = 5)
for (plots in plot_list){
    print(plots)
}
dev.off()