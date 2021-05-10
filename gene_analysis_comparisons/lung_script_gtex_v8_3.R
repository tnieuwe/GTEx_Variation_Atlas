###lung Anaylsis


library(DESeq2)
#setwd("~/Dropbox/projects/gtex/artery")
load(file="gene_analysis_comparisons/data_in/lung-vsd-mean-filtered.rda")
vsdMeanFiltered <- lungVSDMeanFiltered

## calculate variance for each genes
v <- apply(assay(vsdMeanFiltered),1,var)
ind <- which(v > 4)

## heatmap of high variance genes
library(gplots)
#creating highvar subset
vsdHighVar <- assay(vsdMeanFiltered)[ind,]
rownames(vsdHighVar) <- gtabMeanFiltered$gene_name[ind]
colnames(vsdHighVar) <- vsdMeanFiltered$SAMPID
vsdHighVarCentered <- vsdHighVar - rowMeans(vsdHighVar)
# pdf("~/work2/tnieuwe1/data/gtex_v8/diff_samps/lung/lung-high-variance-heatmap.pdf", width=8, height=8)
# hm <- heatmap.2(vsdHighVarCentered, trace="none",col=bluered(20), 
#                 scale="row", margins=c(4,4), cexCol=0.25, cexRow=0.25,
#                 density.info="histogram",key.title="")
#dev.off()

## correlation between high variance genes
library("pheatmap")
library("RColorBrewer")

#Removing duplicated genes, Tim code
dim(vsdHighVar)
vsdHighVar <-  vsdHighVar[!(duplicated(rownames(vsdHighVar))),]
dim(vsdHighVar)

#Gcor is correlations between each genes
gcor <- cor(t(vsdHighVar), method="kendall")

## critical value
#library(pvrank)
library(SuppDists)
ngene <- nrow(vsdHighVar)
nsamp <- ncol(vsdHighVar)
#tauCrit_old <- round(-qrank(0.01/((ngene^2)-ngene)/2, nsamp, index="kendall", approx = "gaussian")$Cq, 2)
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

## heatmap

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
#anno_colors = list("cluster"=c("white", brewer.pal(n-1, "Spectral")))
anno_colors = list("cluster"=c("white", brewer.pal(n-1, "Spectral")))
names(anno_colors$cluster) <- c("", LETTERS[1:(n)-1])
colors <- colorRampPalette( brewer.pal(11, "RdBu")[2:10] )(255)
brks <- seq(0,2,length=length(colors)+1)

# chm <- pheatmap(1-gcor, col=colors, breaks=brks,
#                 cluster_rows=dendo, cluster_col=dendo, 
#                 legend_breaks=c(2,1.5,1,0.5,0), 
#                 legend_labels=c("-1.0","-0.5","0","0.5","1.0"),
#                 main="", fontsize=5, fontsize_row=3, fontsize_col=3,
#                 
#                 annotation_col = anno, annotation_row = anno, 
#                 
#                 annotation_colors = anno_colors,
#                 border_color = NA,
#                 filename="~/work2/tnieuwe1/data/gtex_v8/diff_samps/lung/lung-between-gene-correlation-high-variance-genes.pdf", 
#                 width=10, height=8)
#dev.off()

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



write.csv(gclusts, file="gene_analysis_comparisons/data_out/lung-gene-clusters-high-variance-old.csv",
          quote=FALSE, row.names=FALSE)

##############

## stratify data by inclusion in clusters 
ind <- which(rownames(vsdHighVar) %in% rownames(gcor)[clusts>0])
vsdHighVarClust <- vsdHighVar[ind,]
vsdHighVarCenteredClust <- vsdHighVarCentered[ind,]
vsdHighVarNoise <- vsdHighVar[-ind,]

gcorClust <- gcor[ind,ind]
dendoClust <- hclust(d=as.dist(1-abs(gcorClust)), method = "average")
cnamesClust <- cnames[ind]

## clusters A-H
n <- length(table(cnamesClust))
anno = data.frame("cluster"=factor(cnamesClust))
rownames(anno) = rownames(vsdHighVarClust)
anno_colors = list("cluster"=brewer.pal(n, "Spectral"))
names(anno_colors$cluster) <- LETTERS[1:n]
colors <- colorRampPalette( brewer.pal(11, "RdBu")[2:10] )(255)
brks <- seq(0,2,length=length(colors)+1)
# chm <- pheatmap(1-gcorClust, col=colors, breaks=brks,
#                 cluster_rows=dendoClust, cluster_col=dendoClust, 
#                 legend_breaks=c(2,1.5,1,0.5,0), 
#                 legend_labels=c("-1.0","-0.5","0","0.5","1.0"),
#                 main="", fontsize=5, fontsize_row=3, fontsize_col=3,
#                 annotation_col = anno, annotation_row = anno, 
#                 annotation_colors = anno_colors,
#                 border_color = NA,
#                 filename="~/work2/tnieuwe1/data/gtex_v8/diff_samps/lung/lung-between-gene-correlation-high-variance-genes-clean.pdf", 
#                 width=10, height=8)
#dev.off()

## remake with non A-H genes
gcorNoise <- gcor[-ind,-ind]
dendoNoise <- hclust(d=as.dist(1-abs(gcorNoise)), method = "average")

colors <- colorRampPalette( brewer.pal(11, "RdBu")[2:10] )(255)
brks <- seq(0,2,length=length(colors)+1)
# chm <- pheatmap(1-gcorNoise, col=colors, breaks=brks,
#                 cluster_rows=dendoNoise, cluster_col=dendoNoise, 
#                 legend_breaks=c(2,1.5,1,0.5,0), 
#                 legend_labels=c("-1.0","-0.5","0","0.5","1.0"),
#                 main="", fontsize=10, fontsize_row=10, fontsize_col=10,
#                 border_color = NA,
#                 filename="~/work2/tnieuwe1/data/gtex_v8/diff_samps/lung/lung-between-gene-correlation-high-variance-genes-noise.pdf", 
#                 width=10, height=8)
# #dev.off()


#####################

## add subject specific data
subjtab <- read.table(file="global_in/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", fill=T,
                      stringsAsFactors=FALSE, sep="\t", header=T)
subjIDs <- sapply(strsplit(colnames(vsdHighVarClust),"-"),function(x) paste(x[1],x[2],sep="-"))
map <- match(subjIDs,subjtab$SUBJID)
identical(subjIDs,subjtab$SUBJID[map])

gender <- subjtab$SEX[map]
age <- subjtab$AGE[map]
death <- subjtab$DTHHRDY[map]
death3 <- death
death3[death3>1] <- 2

## add sample specific data
samptab <- read.table("global_in/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                      header=TRUE,stringsAsFactors=FALSE,sep="\t",fill=TRUE,
                      quote="")
map <- match(colnames(vsdHighVarClust),samptab$SAMPID)
identical(colnames(vsdHighVarClust),samptab$SAMPID[map])
samptab <- samptab[map,] 

## extract batch date as a numeric factor
batch <- as.Date(samptab$SMNABTCHD, format="%m/%d/%Y")
batch2 <- as.numeric(as.factor(batch))
bcols <- colorRampPalette( brewer.pal(9, "Purples")[2:9] )(90)

## extract RNA integrity number
rin <- samptab$SMRIN
rcols <- colorRampPalette( brewer.pal(9, "YlGn"))(30)

## heatmap with sample annotation
library(gplots)
library(fields)
library(GMD)
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
#source("heatmap.3.R")
clab <- cbind(rcols[cut(rin, 30)],
              bcols[batch2], 
              c("grey","black")[gender],
              c("dark green","green","light green")[death3+1]
)
colnames(clab) <- c("RIN","Date","Gender","Death")

# pdf("~/work2/tnieuwe1/data/gtex_v8/diff_samps/lung/lung-high-var-heatmap-clean-subj-info.pdf", width=8, height=8)
# hm3 <- heatmap.3(vsdHighVarCenteredClust[order(cnamesClust, decreasing = TRUE),], 
#                  Rowv=NA, trace="none",
#                  col=bluered(20),
#                  scale="none", 
#                  margins=c(2,12), labCol=NA, labRow=NA,
#                  density.info="none", KeyValueName="Residuals", keysize=1,
#                  xlab="tish Samples", ylab="", dendrogram="none",
#                  sidesize=4, 
#                  ColSideColors=clab, ColSideColorsSize=4 
#                  #add.exprs=mtext(side=2,paste("Cluster", LETTERS[1:8]), at=c(0.15,0.47,0.68,0.8,0.88,0.95,1.01),line=2, las=2)
# )
##dev.off()
# 
# legend(x=0.83,y=0.9,c("Ventilator","Violent","Natural"), title="Death",
#        fill=c("dark green","green","light green"), border=FALSE, bty="n", cex=1)
# legend(x=0.83,y=0.7,c("Female","Male"), title="Gender",
#        fill=c("black","grey"), border=FALSE, bty="n", cex=1)
# legend(x=0.83,y=0.55,c("06/27/11","","09/12/13","","01/28/14"), title="Date",
#        fill=rep("white",5), border=FALSE, bty="n", cex=1)
# colorbar.plot(x=0.87, y=0.42, strip=sort(unique(batch2), decreasing=TRUE), col=bcols, 
#               strip.width=0.02, strip.length=0.14, horizontal=FALSE)
# legend(x=0.83,y=0.31,c("3.2","","7.0","","9.4"), title="  RIN",
#        fill=rep("white",5), border=FALSE, bty="n", cex=1)
# colorbar.plot(x=0.87, y=0.18, strip=seq(from=max(rin), to=min(rin), length=length(rin)), 
#               col=rcols, strip.width=0.02, strip.length=0.14, horizontal=FALSE)
#dev.off()

save(vsdHighVar, vsdHighVarCentered, vsdHighVarClust, vsdHighVarNoise, 
     vsdHighVarCenteredClust, cnames, cnamesClust,
     file="~/work2/tnieuwe1/data/gtex_v8/diff_samps/lung/lung-high-variance-gene-data.rda")

#########################################

## compute avg cluster expression for each sample
#load(file="~/work2/tnieuwe1/data/gtex_v8/diff_samps/lung/lung-high-variance-gene-data.rda")

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

save(clusterProfiles, file="~/work2/tnieuwe1/data/gtex_v8/diff_samps/lung/lung-cluster-profiles.rda")
write.csv(clusterProfiles, file="~/work2/tnieuwe1/data/gtex_v8/diff_samps/lung/lung-cluster-profiles.csv", quote=FALSE)
