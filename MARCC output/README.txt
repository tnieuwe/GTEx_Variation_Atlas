Folders:
*old_files*: Contains older versions of MARCC pipeline outputs
either because we changed parameters or the pipeline itself
changed.
*variance_genes_(spearman|kendall)*: For each tissue the genes
that clustered in their respective clusters, for either a spearman
or kendall correlation.
*variance_heatmaps_(spearman|kendall)*: For each tissue the
corresponding pheatmap of their respective correlation matrix.
*variance_profiles_(spearman|kendall)*: For each tissue, the
calculated Z-score association between each samples and cluster.
*all_tissues_umaps*: A file including PDFs of UMAPs from each
tissue. The UMAPs are derived from the high variance genes
and they are colored on the Z-scores of various clusters.
For the most part contamination and sex clusters are not included
in this analysis. 
Changes per date:
2021-5-19:
+ mcall_cluster_generator
	+ sub_clust_off_cut = 0.8 from 0.85
	+ agglo_cut_off	= 0.7 from 0.85
2021-6-14:
+ Umap seed set to 123