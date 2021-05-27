
Last updated 2021-5-27
FILES:
across_cluster_gene_correlation.rmd:
This read in the cluster gene results from the MARCC pipeline
into a a nested list with tissues at the top layer,
followed by clusters, and finally the genes in the cluster.
We then filter this data down to gene that appear >6 times the
dataset as high variance. From there we then "correlate" the
genes by testing how often they show up together in clusters.

Using this "pseudo correlation" we then cluster the values
using euclidian distance. From there, by hand, we determine
what seem to be clusters of interest. We refer to these as
meta-clusters.

With these metaclusters we again go through all the clusters for
all tissues and determine what percentage of the clusters are
explained by meta clusters. From there create a plot for each
tissue, saved under the file all_tissue_cluster_composition.pdf.

The next step is to create a write out that will be used for
further steps in automatting pipeline discovery. This is still
to be discussed on the optimal method to do so.
