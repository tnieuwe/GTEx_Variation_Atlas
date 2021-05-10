load("~/work2/tnieuwe1/data/gtex_v7/gtex-gene-counts-v7.rda")

## analyze GTEx general data to look for variability between samples
#setwd("~/tnieuwe1/projects/gtex/general")
load(file="~/work2/tnieuwe1/data/gtex_v7/diff_samps/general/gtex-gene-counts-general.rda")

## transform data to a DESeq object
library(DESeq2)

###This should be obsolete as they're alread split
## split based on coronary or tibial general
#iTibial <- which(generaltab$SMTSD == "Artery - Tibial")

#Original code doesn't work because colnames are the same, below is my fix
rownames(generaltab) <- generaltab$SAMPID

generalDES <- DESeqDataSetFromMatrix(countData = generaldat,
                                            colData = generaltab,
                                            design = ~ 1)
rownames(generalDES) <- gtab$Name


## variance stabilizing transformation
generalSD <- varianceStabilizingTransformation(generalDES)
save(generalSD, gtab, file="~/work2/tnieuwe1/data/gtex_v7/diff_samps/general/general-vsd-unfiltered.rda")

##############################################################

## filter low expressed genes - general
m <- rowMeans(assay(generalSD))
pdf("~/work2/tnieuwe1/data/gtex_v7/diff_samps/general/general-avg-vsd-histogram.pdf", width=8, height=8)
hist(m, breaks=25, main="", xlab="Normalized transformed counts")
abline(v=5, lwd=3, lty=2)
dev.off()

## select clearly expressed genes
ind <- which(m>5)
gtabMeanFiltered <- gtab[ind,]
generalVSDMeanFiltered <- generalSD[ind,]
save(gtabMeanFiltered, generalVSDMeanFiltered, file="~/work2/tnieuwe1/data/gtex_v7/diff_samps/general/general-vsd-mean-filtered.rda")

## filter raw count DESeq object
generalDESMeanFiltered <- generalDES[ind,]
identical(rownames(generalDESMeanFiltered), gtabMeanFiltered$Name)
save(generalDESMeanFiltered, gtabMeanFiltered, file="~/work2/tnieuwe1/data/gtex_v7/diff_samps/general/general-counts-mean-filtered.rda")

##############################################################

