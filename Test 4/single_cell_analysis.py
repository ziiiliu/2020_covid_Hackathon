import pandas as pd
import numpy as np
import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, facecolor='white')

df_iter = pd.read_csv('covid_scRNA_data.csv', chunksize = 1000, iterator = True)
anno = pd.read_csv('covid_scRNA_metadata.csv')

for iter_num, chunk in enumerate(df_iter, 1):
    print('Processing iteration' + str(iter_num))
    adata = chunk
    adata.obs['cluster'] = anno['cell_cluster']
    adata.obs['cluster'] = adata.obs['cluster'].astype('category')

#sc.pl.tsne(adata, color = ['cluster'])

sc.pl.higeset_expr_genes(adata, n_top=20)
#sc.pp.filter_cells(adata, min_genes =200)
#sc.pp.filter_genes(adata, min_cells =3)

sc.pl.umap(adata,color = ['leiden', 'CST3', 'NKG7'])
adata.write(results_file)
sc.tl.rank_genes_groups(adata, 'leiden', method = 't-test)
sc.pl.rank_genes_groups(adata, n_genes = 25, sharey = False)
