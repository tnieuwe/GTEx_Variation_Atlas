
#general version

#load datasets (heavy part of code)
load("~/work2/tnieuwe1/data/gtex_v8/gtex-gene-counts-v8.rda")


## select desired data
#ind <- which(stab$SMTSD%in%c("tish"))
ind <- which(stab$SMTSD%in%c("tish"))

generaldat <- dat[,ind]
generaltab <- stab[ind,]

## load subject annotation
subjtab <- read.table("~/work2/tnieuwe1/data/gtex_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                      header=TRUE,stringsAsFactors=FALSE,sep="\t",fill=TRUE,
                      quote="")

## match subject annotation to sample annotation
subjIDs <- sapply(strsplit(generaltab$SAMPID,"-"),
                  function(x) paste(x[1],x[2],sep="-"))
map <- match(subjIDs,subjtab$SUBJID)
identical(subjIDs,subjtab$SUBJID[map])

## merge subject and sample annotation
generaltab <- data.frame(subjtab[map,], generaltab)

save(gtab, generaldat, generaltab, file="~/work2/tnieuwe1/data/gtex_v8/diff_samps/general/gtex-gene-counts-general.rda")
write.csv(generaltab, file="~/work2/tnieuwe1/projects/gtex_v8/general/gtex_general.csv", row.names=FALSE, quote=TRUE)
