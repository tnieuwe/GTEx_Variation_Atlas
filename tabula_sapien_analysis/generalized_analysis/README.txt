## Purpose
The purpose of this code is to generate for each tissue relevant cell markers for each cluster

## Plan
1. Download the figshare files from their respective location into:
	work2/tniewe1/data/gtex_v8/tabula_sapiens/ (downloaded 29/07/2021)
"wget https://ndownloader.figshare.com/articles/14267219/versions/3" (ran 29/7/21)
2. For each tabula organ pull all gtex tissues that are similair to it.
3. Pull out the gene names from the gtex clusters to focus the analysis on a smaller subset
	of genes.
4. Run the ROC method and create for each organ its own output dataframe of marker genes
	labelled by cell type by also what tabula tissue its coming from
5. Run the normal method too to generate associated p-values, join this will ROC method for
	the purpose of removing non-sig findings as a very low level filter step.
5. Collate all of this data into a master list to generate and eventual equation with.