
#lung version

#load datasets (heavy part of code)
load("gene_analysis_comparisons/original_lung_analysis/data_in/gtex-gene-counts-v8-old.rda")


## select desired data
#ind <- which(stab$SMTSD%in%c("Lung"))
ind <- which(stab$SMTSD%in%c("Lung"))

lungdat <- dat[,ind]
lungtab <- stab[ind,]

## load subject annotation
subjtab <- read.table("global_in//GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                      header=TRUE,stringsAsFactors=FALSE,sep="\t",fill=TRUE,
                      quote="")

## match subject annotation to sample annotation
subjIDs <- sapply(strsplit(lungtab$SAMPID,"-"),
                  function(x) paste(x[1],x[2],sep="-"))
map <- match(subjIDs,subjtab$SUBJID)
identical(subjIDs,subjtab$SUBJID[map])

## merge subject and sample annotation
lungtab <- data.frame(subjtab[map,], lungtab)

save(gtab, lungdat, lungtab, file="gene_analysis_comparisons/original_lung_analysis/data_out/gtex-gene-counts-lung.rda")
write.csv(lungtab, file="gene_analysis_comparisons/original_lung_analysis/data_out/gtex_lung.csv", row.names=FALSE, quote=TRUE)
