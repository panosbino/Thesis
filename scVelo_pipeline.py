import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

sample_one = anndata.read_loom("subset_final_possorted_bam_HF5H3.loom")
sample_obs = pd.read_csv("cellID_obs.csv")
umap_cord = pd.read_csv("cell_embeddings.csv")
cell_clusters = pd.read_csv("clusters.csv")
#umap = pd.read_csv("umap.csv")

my_list= []

for i in range(0,1289):
    my_list.append(sample_one.obs.index.str.split(':')[i][1]+"-1")

sample_one.obs.index = my_list

sample_one.obs.index
#Let's cast our index as a data frame and change the column name

sample_one_index = pd.DataFrame(sample_one.obs.index)
sample_one_index = sample_one_index.rename(columns = {0:'Cell ID'})

umap = umap_cord.rename(columns = {'Unnamed: 0':'Cell ID'})

umap_ordered = sample_one_index.merge(umap, on = "Cell ID")

umap_ordered = umap_ordered.iloc[:,1:]
sample_one.obsm['X_umap'] = umap_ordered.values
cell_clusters = cell_clusters.rename(columns = {'rowname':'Cell ID'})
clusters_ordered = sample_one_index.merge(cell_clusters, on = "Cell ID")
clusters = clusters_ordered.iloc[:,2:]
sample_one.uns['Cluster_colors'] = clusters.values
scv.pp.filter_and_normalize(sample_one)
scv.pp.moments(sample_one)
scv.tl.velocity(sample_one, mode = "stochastic")
scv.tl.velocity_graph(sample_one)
scv.pl.velocity_embedding(sample_one, basis = 'umap')

scv.pl.velocity_embedding_stream(sample_one, basis='umap', color = sample_one.uns['Cluster_colors'], save = "velo_stream.svg")
scv.pl.proportions(sample_one)
scv.pl.velocity_embedding(sample_one, arrow_length=3, arrow_size=2, dpi=120, color = sample_one.uns['Cluster_colors'],save = "velo_arrows.svg")
