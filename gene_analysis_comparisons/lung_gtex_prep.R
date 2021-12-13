## This file prepares all lung samples to be analyzed by later pipelienes
## this section is priamrilt he subsetting out of lung data and its normalization
## using DESeq2's VST transformation.

##Input variables you care about
cur_var <- 2.5
numeric_var <- cur_var
# if (cur_var %% 2 != 0) {
#     cur_var <- paste0(cur_var, "n-half")
# }
#cur_var <- 


## Load in libraries
library(DESeq2)
library(tidyverse)
library(SuppDists)


load("global_in/gtex-gene-counts-v8.rda")
#load("gene_analysis_comparisons/data_in/gtex-gene-counts-v8-old.rda")
#version <- "pre_fix"
version <- ""
tissue <- "Artery - Coronary"
source("global_in/general_scripts.R")
pheno_dat <- read.csv(file = "global_in/gtex_phenotypes_v8.csv")

## Filter down to lung data

## Check if true for indexing
(sum(stab$SAMPID == colnames(dat))) == ncol(dat)
## Filter
samp_ind <- (stab$SMTSD == tissue)
dat_lung <- dat[,samp_ind]
stab_lung <- stab[samp_ind,]

## Make DESeq object
des_dat <- DESeqDataSetFromMatrix(countData = dat_lung,
                                  colData = stab_lung,
                                  design = ~1)
## Normalize data
#vst_lung <- vst(des_dat)
vst_lung <- varianceStabilizingTransformation(des_dat)

# system.time(vst(des_dat))
# system.time(varianceStabilizingTransformation(des_dat))

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
ind <- which(v > numeric_var)
high_var_lung <- vst_lung_mean_filtered[ind,]
high_var_gtab <- gtab_lung[ind,]
high_var_matrix_lung <- assay(high_var_lung)

save(vst_lung_mean_filtered, stab_lung, gtab_lung,
     high_var_lung,
     high_var_gtab,
     file = paste0("gene_analysis_comparisons/data_out/lung_vst_var_",cur_var,".rda"))
save(vst_lung_mean_filtered, stab_lung, gtab_lung,
     high_var_lung,
     high_var_gtab,
     file =  paste0("gene_analysis_comparisons/data_in/lung_vst_var_",cur_var,".rda"))

## Alphnum clustering method
rownames(high_var_matrix_lung) <- high_var_gtab$gene_name
vsd_high_var_centered_lung <- high_var_matrix_lung - rowMeans(high_var_matrix_lung)

gcor <- cor(t(high_var_matrix_lung))

## Generate critical value
ngene <- nrow(high_var_matrix_lung)
nsamp <- ncol(high_var_matrix_lung)
tauCrit <- round(qKendall(0.01/((ngene^2)-ngene)/2, nsamp, lower.tail = F), 2)

## cluster based on absolute correlation distance
dendo <- hclust(d=as.dist(1-abs(gcor)), method = "average")

## clean up heatmap to focus on interesting gene clusters
clusts <- cutree(dendo, h=1-tauCrit)
tabclusts <- table(clusts)
ind <- which(tabclusts >= 6) #arbitrary
clusts[!(clusts %in% names(tabclusts)[ind])] <- 0
ntmp <- names(clusts)
clusts <- as.numeric(as.factor(clusts))-1
names(clusts) <- ntmp


#Naming clusters based on their clust values
#Very manual process
n <- length(table(clusts))
cnames <- rep("", length(clusts))

#Tim for loop
for (clst in 1:(n)){
    cnames[clusts==clst] <- LETTERS[clst]
}

anno = data.frame("cluster"=factor(cnames))
rownames(anno) = names(clusts)
## make matrix to output
gclusts <- matrix("", nrow=max(table(cnames)), ncol=length(table(cnames)))
colnames(gclusts) <- levels(anno$cluster)
for(k in 1:length(table(cnames))){
    tmp <- colnames(gclusts)[k]
    gclusts[1:sum(cnames==tmp), k] <- rownames(gcor)[cnames==tmp]
}

#Sub clusters code! Now annotated

for (column in LETTERS[1:(length(unique(cnames)) -1 )]) {
    #Pull out the first column
    test_col <- gclusts[,column]
    first_gene <- test_col[1]
    #Make an index to subset gcor to just a single clusters
    ind <- which(rownames(gcor) %in% test_col)
    sub_cor <- gcor[ind,ind]
    #Test if the cluster has any negatively correlating genes
    ind2 <- sub_cor[,first_gene] < 0
    
    #If there are negatively correlating genes...
    if (sum(ind2) > 0) {
        #generate sub cluster 1
        sub_1 <- rownames(sub_cor)[!ind2]
        #Subcluster 2
        sub_2 <- rownames(sub_cor)[ind2]
        
        #Remover old column to replace with new columns
        gclusts <- as.matrix(gclusts[,!(colnames(gclusts) == column)])
        
        
        #nrow(gclusts) Not sure if useful, keeping for now
        #add sub-columns
        
        #Make a single matrix of same rows as gclusts, one for each sub-cluster. 
        matrix_1  <- as.matrix(c(sub_1 ,rep("", nrow(gclusts) - length(sub_1))), )
        colnames(matrix_1) <- paste0(column,"-1")
        
        # cbind(gclusts, matrix_1, matrix_2)
        
        
        matrix_2  <- as.matrix(c(sub_2 ,rep("", nrow(gclusts) - length(sub_2))), )
        colnames(matrix_2) <- paste0(column,"-2")
        
        #Add subclusters to gclusts
        gclusts <-   cbind(gclusts, matrix_1, matrix_2)
        
    }
    
    
}



#Realphebatize columns for sake of beauty 
gclusts <- gclusts[,order(colnames(gclusts))]

write.csv(gclusts, file=paste0("gene_analysis_comparisons/data_out/lung-gene-clusters-",cur_var,"-variance.csv"),
          quote=FALSE, row.names=FALSE)



### Generate profiles

## extract expression for each cluster, convert to zscores, and summarize profile
clusterProfiles <- matrix(nrow=ncol(vsd_high_var_centered_lung), ncol=(ncol(gclusts) - 1))
for(k in 1:ncol(clusterProfiles)){
    # ind <- which(cnamesClust == LETTERS[k])
    
    #New code picks from gclusts
    ind <- rownames(vsd_high_var_centered_lung) %in%  gclusts[,1 + k]
    
    #FOr normal clusters 
    if (length(unique(gclusts[,1 + k])) != 2) {
        
        etmp <- vsd_high_var_centered_lung[ind,]
        ztmp <- etmp / apply(etmp, 1, mad)
        #Below removes inf values
        ztmp <- ztmp[is.finite(rowSums(ztmp)),]
        
        clusterProfiles[,k] <- colMeans(ztmp)
    } else {
        #For single genes such as XIST
        etmp <- vsd_high_var_centered_lung[ind,]
        clusterProfiles[,k] <- etmp
        
        
    }
    
    
    
    
}
rownames(clusterProfiles) <- colnames(vsd_high_var_centered_lung)
colnames(clusterProfiles) <- colnames((gclusts))[-1]

#save(clusterProfiles, file="gene_analysis_comparisons/original_lung_analysis/data_out/lung-cluster-profiles.rda")
write.csv(clusterProfiles,
          file= paste0("gene_analysis_comparisons/data_out/lung-cluster-profiles-",cur_var,"-var.csv"), quote=FALSE)
