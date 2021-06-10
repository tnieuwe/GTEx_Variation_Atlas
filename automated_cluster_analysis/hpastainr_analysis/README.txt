files:
*hpastainr_analysis.rmd*: An rmarkdown file that runs all
clusters through the HPAStainR R package to dtermine if clusters
are associated with specific cell types.
Last ran: 6/10/2021
hpar_1.34.0
HPAStainR_1.3.1

folders: 
**data_in**:
+ *ordered_all_tiss_hpa_names.csv*: A file connecting HPA organs to
GTEx organs. Curated by Tim N.
**data_out**: 
+ *hpastainr_analysis_uncleaned.rda*: An rda file containing
all results from the running of HPAStainR on the clusters in an
R list format. This is an intermediate file saved due to the
long run time of HPAStainR.
+ *<tissue>_hpastainr_results.csv*: The filtered results from
the HPA analysis including only data from matched GTEx HPA
tissues. 
