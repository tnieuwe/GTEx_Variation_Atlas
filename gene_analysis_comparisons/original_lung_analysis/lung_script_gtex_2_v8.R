#load("~/work2/tnieuwe1/data/gtex_v8/gtex-gene-counts-v8.rda")

## analyze GTEx lung data to look for variability between samples
#setwd("~/tnieuwe1/projects/gtex/lung")
load(file="gene_analysis_comparisons/original_lung_analysis/data_out/gtex-gene-counts-lung.rda")

## transform data to a DESeq object
library(DESeq2)

###This should be obsolete as they're alread split
## split based on coronary or tibial lung
#iTibial <- which(lungtab$SMTSD == "Artery - Tibial")

#Original code doesn't work because colnames are the same, below is my fix
rownames(lungtab) <- lungtab$SAMPID

lungDES <- DESeqDataSetFromMatrix(countData = lungdat,
                                            colData = lungtab,
                                            design = ~ 1)
rownames(lungDES) <- gtab$Name


## variance stabilizing transformation
lungVSD <- varianceStabilizingTransformation(lungDES)
save(lungVSD, gtab, file="gene_analysis_comparisons/original_lung_analysis/data_out/lung-vsd-unfiltered.rda")

##############################################################

## filter low expressed genes - lung
m <- rowMeans(assay(lungVSD))
pdf("gene_analysis_comparisons/original_lung_analysis/data_out/lung-avg-vsd-histogram.pdf", width=8, height=8)
hist(m, breaks=25, main="", xlab="Normalized transformed counts")
abline(v=5, lwd=3, lty=2)
dev.off()

## select clearly expressed genes
ind <- which(m>5)
gtabMeanFiltered <- gtab[ind,]
lungVSDMeanFiltered <- lungVSD[ind,]
save(gtabMeanFiltered, lungVSDMeanFiltered, file="gene_analysis_comparisons/original_lung_analysis/data_out/lung-vsd-mean-filtered.rda")

## filter raw count DESeq object
lungDESMeanFiltered <- lungDES[ind,]
identical(rownames(lungDESMeanFiltered), gtabMeanFiltered$Name)
save(lungDESMeanFiltered, gtabMeanFiltered, file="gene_analysis_comparisons/original_lung_analysis/data_out/lung-counts-mean-filtered.rda")

##############################################################

