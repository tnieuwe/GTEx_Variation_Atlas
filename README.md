# GTEx_Variation_Atlas

## Overall goal of project
Determining not only what variation exists in GTEx, but what drives said variation. 


## Data acquistion:
+ gtex_read_counts.gct* The version 8 readcounts, are downloaded through R scripts, for most recent download
see gtex_read_counts.gct last date modified in Setup Code/downloaded_input/gtex_read_counts.gct
+ *gtex_phenotypes_v8*: A file containing deeper phenotype data on GTEx participants
was originally acquired through dbGAP on 2019-10-21
+ *GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt*: a file containing sample
data on GTEx samples was last download 2019-9-30
+ GTEx_Analysis_v8_Annotations_SubjectAttributesDS.txt*: A file containing subject
data on GTEx subjects was last downloaded 2019-8-26.


## File descriptiojns

*gen_script_gtex_v8_1*:
1. Separates tissues readcount matrices from one another.

*gen_script_gtex_v8_2*:

2. Normalizes the data using DESEq2's VST() function before filtering genes on a
threshold that rowMeans() should be >5 normalized counts.

*gen_script_gtex_v8_3_v2*:

**mccall_analysis_corr_step.R**:

3. The gene data is further filtered to only genes in the .98 percentile of high
variance genes for each tissue.
4. Then using a Kendall's tau correlation the high variance genes generate a correlation matrix.

**mccall_cluster_generator.R**:

5. This correlation matrix is used to create euclidean distances to cluster the data on absolute values.
Generate actual clusters through a kendall's tau critical value.
6. Further break up clusters based on how tightly they correlate using a quantile
threshold of correlation, and new agglometric clustering method.

**sub_cluster_maker.R**:

7. Break the clusters up on if they positively or negatively correlate with one another.

**mccall_profile_maker.R**:

8. Using the genes from each cluster, generate a Z-score of samples association with 
each cluster using the normalized count data.
