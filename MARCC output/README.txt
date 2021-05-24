Folders:
*old_files*: Contains older versions of MARCC pipeline outputs
either because we changed parameters or the pipeline itself
changed.
*variance_genes_(spearman|kendall)*: For each tissue the genes
that clustered in their respective clusters, for either a spearman
or kendall correlation.
*variance_heatmaps_(spearman|kendall)*: For each tissue the
corresponding pheatmap of their respective correlation matrix.
*variance_profiles_(spearman|kendall)*: FOr each tissue, the
calculated Z-score association between each samples and cluster.
Changes per date:
2021-5-19:
+ mcall_cluster_generator
	+ sub_clust_off_cut = 0.8 from 0.85
	+ agglo_cut_off	= 0.7 from 0.85