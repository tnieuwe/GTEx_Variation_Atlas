## This file prepares all lung samples to be analyzed by later pipelienes
## this section is priamrilt he subsetting out of lung data and its normalization
## using DESeq2's VST transformation.

## Load in libraries
library(DESeq2)
library(tidyverse)


load("global_in/gtex-gene-counts-v8.rda")
source("global_in/general_scripts.R")
pheno_dat <- read.csv(file = "global_in/gtex_phenotypes_v8.csv")

## Filter down to lung data

## Check if true for indexing
(sum(stab$SAMPID == colnames(dat))) == ncol(dat)
## Filter
samp_ind <- (stab$SMTSD == "Lung")
dat_lung <- dat[,samp_ind]
stab_lung <- stab[samp_ind,]

## Make DESeq object
des_dat <- DESeqDataSetFromMatrix(countData = dat_lung,
                                  colData = stab_lung,
                                  design = ~1)
## Normalize data
vst_lung <- vst(des_dat)

## filter low expressed genes 
m <- rowMeans(assay(vst_lung))
hist(m, breaks=25, main="", xlab="Normalized transformed counts")
abline(v=5, lwd=3, lty=2)
dev.off()

## select clearly expressed genes
ind <- which(m>5)
gtab_lung <- gtab[ind,]
vst_lung_mean_filtered <- vst_lung[ind,]

## Join stab with pheno to combine all phenotype data
stab_lung$SUBJID <- GTEx_SAMPID_to_SUBJIUD(stab_lung$SAMPID)

stab_lung <- left_join(stab_lung, pheno_dat, by = "SUBJID") %>%
    select(SAMPID, SUBJID, everything())

# save(vst_lung_mean_filtered, stab_lung, gtab_lung, file = "gene_analysis_comparisons/data_out/lung_vst.rda")
# save(vst_lung_mean_filtered, stab_lung, gtab_lung, file = "gene_analysis_comparisons/data_in/lung_vst.rda")

## Filter down to high variance genes
v <- apply(assay(vst_lung_mean_filtered),1,var)
ind <- which(v > 4)
high_var_lung <- vst_lung_mean_filtered[ind,]
high_var_gtab <- gtab_lung[ind,]

save(vst_lung_mean_filtered, stab_lung, gtab_lung,
     high_var_lung,
     high_var_gtab,
     file = "gene_analysis_comparisons/data_out/lung_vst.rda")
save(vst_lung_mean_filtered, stab_lung, gtab_lung,
     high_var_lung,
     high_var_gtab,
     file = "gene_analysis_comparisons/data_in/lung_vst.rda")
