### general analysis

## Load libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(SuppDists)
library(reshape2)

## Load function
source("~/work2/tnieuwe1/data/gtex_v8/pipeline_scripts/mccall_analysis_corr_step.R")
source("~/work2/tnieuwe1/data/gtex_v8/pipeline_scripts/mccall_cluster_generator.R")
source("~/work2/tnieuwe1/data/gtex_v8/pipeline_scripts/mccall_profile_maker.R")
source("~/work2/tnieuwe1/data/gtex_v8/pipeline_scripts/sub_cluster_maker.R")
source("~/work2/tnieuwe1/data/gtex_v8/pipeline_scripts/sub_cluster_maker.R")
source("~/work2/tnieuwe1/data/gtex_v8/pipeline_scripts/unfiltered_heatmap.R")

## Set up strings
profiles_path <- "~/work2/tnieuwe1/data/gtex_v8/diff_samps/general/general-cluster-profiles.csv"
gene_clust_path <- "~/work2/tnieuwe1/data/gtex_v8/diff_samps/general/general-gene-clusters-high-variance.csv"
heatmap_path <- paste0("~/work2/tnieuwe1/data/gtex_v8/",
                       "diff_samps/general/general-between-gene",
                       "-correlation-high-variance-genes.pdf")


## Load data
# load(file="gene_analysis_comparisons/data_in/general-vsd-mean-filtered.rda")
# vsdMeanFiltered <- generalVSDMeanFiltered
pipeline_dat <- mccall_analysis_corr_step("general",
                                     path_to_files = "~/work2/tnieuwe1/data/gtex_v8/diff_samps/general/",
                                     variance_type = "quantile",
                                     variance_value = .98,
                                     correlation = "spearman",
                                     MARCC = TRUE)
pipeline_dat <- mccall_cluster_generator(pipeline_dat,
                                         cluster_breaker = TRUE)
pipeline_dat <- sub_cluster_maker(pipeline_dat)
pipeline_dat <- mccall_profile_maker(pipeline_dat)

## Make figures and save data
unfiltered_heatmap(pipeline_dat,
                   save_file = heatmap_path) 
write.csv(pipeline_dat[["sample_profiles"]], file=profiles_path, quote=FALSE)
write.csv(pipeline_dat[["cluster_genes"]], file = gene_clust_path)