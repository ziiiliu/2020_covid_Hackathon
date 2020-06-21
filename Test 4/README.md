# Single-cell RNA expression from infected patients

This dataset measures the transcriptomes of single-cells isolated from the fluid milieu of
the lungs (BALF). This lung fluid contains a variety of immune cells that circulate through 
the lungs as well as epithelial cells that make up the lungs which can be carried
away in the fluid as well. Cells are derived from 12 different patients : 6 healthy controls,
3 patients with moderate COVID-19, and 6 patients with severe COVID-19. 

To help make the data easier to analyze, a coarse-grained analysis has already been done
to assign what cell type each single cell most closely matches. In the metadata, the column
for "cell_cluster" groups the cells into 30 clusters of closely-related cells, as determined 
by unsupervised clustering algorithms (leiden). The approximate cell types that each cluster
may correspond to are indicated in the "Celltype_tentative" column. Note some celltypes made be 
made of multiple clusters if there is some element that's distinct in what genes those 
subclusters are expressing. 

Additional quality-control metrics are also pre-calculated for you for each cell, including the total
number of distinct genes/transcripts detected per cell, total number of RNA molecules detected, and 
percent of all RNAs detected that derive from the mitochondria. All cells have already passed a
standard quality-control pipeline. Expression values in the matrix have already been normalized and
log-transformed according to a typical data-processing pipeline. 