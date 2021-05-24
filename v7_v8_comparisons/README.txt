FOLDERS:
*input*: any files required for input, currently empty.
*output*: The result of the analyses from cluster_version_comparison.rmd
*run_history*: Includes the previous full runs of the analysis
and theirsession infos.
FILES:
*cluster_version_comparison.rmd*: Compares the oriignal outputs
of GTEx v7 and GTEx v8, note the data used here was generated 
before the creation of the new McCall pipeline and lacks the
agglo-cluster approach and reactive thresholds. 

FILES IN *output*:
*cluster_n_version_comparison*: A barchart comparing the
differences in how many clusters were recorded in the pipeline
on GTEx V7 compared to GTEx V8.
*gene_unique_to_analysis*: A barchart comparing how many genes
are unique to either V7 or V8 of GTEx.
*gene_unique_to_analysis_no_X*: Same plot as above,
but non-clustering genes (cluster X) is removed.
