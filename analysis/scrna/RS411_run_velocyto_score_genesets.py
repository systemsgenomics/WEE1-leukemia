#RS4;11 scRNA-seq velocity analysis with dynamical model and gene set scoring

#Running RNA velocity analysis as presented in scvelo package documentation 
#Sample S6 = treatment, S8 = ctrl
#Velocity is run using dynamical model instead of deterministic

from IPython.core.display import display, HTML
display(HTML("<style>.container { width:90% !important; }</style>"))

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
from matplotlib.colors import ListedColormap
import seaborn as sb

import re
import scvelo as scv
from scanpy import read_10x_h5
import scanpy as sc
import anndata as ad

#scv.logging.print_versions()
#sc.logging.print_versions()
## Running scvelo 0.2.1 (python 3.8.3) on 2020-08-28 17:40.
## scanpy==1.5.1 anndata==0.7.3 umap==0.4.6 numpy==1.18.5 scipy==1.5.0 pandas==1.0.5 scikit-learn==0.23.1 statsmodels==0.11.1 python-igraph==0.8.2 louvain==0.6.1 leidenalg==0.8.1

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, dpi_save=300, vector_friendly=True)

# RS411_velocity_analysis-dynamic-new.h5ad

# A nice color scheme for visualizing gene expression
colors_0 = plt.cm.OrRd(np.linspace(0.05, 1, 128))
colors_1 = plt.cm.Greys_r(np.linspace(0.8,0.9,20))
colors_Comb = np.vstack([colors_1, colors_0])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colors_Comb)

colors_2 = plt.cm.Blues(np.linspace(0.7, 0.05, 128))
colors_3 = plt.cm.Greys_r(np.linspace(0.9,0.8,20))
colors_Comb1 = np.vstack([colors_2, colors_3])
mymap1 = colors.LinearSegmentedColormap.from_list('my_colormap', colors_Comb1)

## INPUT DATA
# Read count data
S6 = read_10x_h5('cellranger/S6/outs/filtered_feature_bc_matrix.h5')
S8 = read_10x_h5('cellranger/S8/outs/filtered_feature_bc_matrix.h5')

# Read spliced/unspliced count data
S6s = scv.read('cellranger/S6/coding_velocyto_res/S6.loom')
S8s = scv.read('cellranger/S8/coding_velocyto_res/S8.loom')

# Get wanted barcodes from text file to filter data
file = open('cellbarcodes_pt_times_RSfilt.txt','r')
barcodes = pd.read_csv('cellbarcodes_pt_times_RSfilt.txt', sep="\t", header=None)
file.close()

# Change suffix to match the adata batches
barcodes[0] = [re.sub('-1_2','_S6',x) for x in barcodes[0]]
barcodes[0] = [re.sub('-1_1','_S8',x) for x in barcodes[0]]

pt_times = barcodes[1]
pt_times_list = list(pt_times)

# Edit cellIDs to match them
S6.obs_names = [re.sub('-1', '', x) for x in S6.obs_names]
S6s.obs_names = [re.sub('.*:|x', '', x) for x in S6s.obs_names]
S8.obs_names = [re.sub('-1', '', x) for x in S8.obs_names]
S8s.obs_names = [re.sub('.*:|x', '', x) for x in S8s.obs_names]

# Add spliced/unspliced counts as layers
S6 = scv.utils.merge(S6, S6s)
S8 = scv.utils.merge(S8, S8s)

# Combine datas
adata = S6.concatenate(S8, batch_categories=['S6', 'S8'], index_unique='_')

# Clean up AnnData object
scv.utils.cleanup(adata)
del adata.obs['Clusters'], adata.obs['_X'], adata.obs['_Y'], adata.obs['initial_size']

# Take subset of adata including only cells that are wanted to the analysis
adata = adata[np.isin(adata.obs_names,barcodes[0])]

# Remove mouse genes
adata = adata[:,adata.var_names.str.startswith('hg19')]
adata.var_names = [re.sub('hg19_', '', x) for x in adata.var_names]

# Save raw counts
adata.raw = scv.pp.log1p(adata, copy=True)

scv.pl.proportions(adata)

# Prelim filter
sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=10)

mito_genes = [name for name in adata.var_names if name.startswith('MT-')]
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
ribo_genes = [name for name in adata.var_names if name.startswith('RP')]
adata.obs['percent_ribo'] = np.sum(adata[:, ribo_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

# Check QC metrics
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito', 'percent_ribo'], stripplot=False, multi_panel=True, groupby='batch')

#From scvelo tutorial 'RNA Velocity Basics':
#All of this is summarized in a single function scv.pp.filter_and_normalize, which essentially runs the following: scv.pp.filter_genes(adata, min_shared_counts=20)
#scv.pp.normalize_per_cell(adata)
#scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
#scv.pp.log1p(adata) Further, we need the first and second order moments (means and uncentered variances) computed among nearest neighbors in PCA space, summarized in ```scv.pp.moments```, which internally computes `scv.pp.pca` and `scv.pp.neighbors`. First order is needed for deterministic velocity estimation, while stochastic estimation also requires second order moments.

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# Read cell cycle genes from gene list by Regev lab
s_genes = pd.read_csv('regev_lab_cell_cycle_genes_S_phase.txt', header=None)[0]
g2m_genes = pd.read_csv('regev_lab_cell_cycle_genes_G2M_phase.txt', header=None)[0]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

#Separate by treatment for downstream analysis
S6 = adata[adata.obs.batch == 'S6']
S8 = adata[adata.obs.batch == 'S8']

sc.pp.pca(S6, n_comps=20, svd_solver='randomized') # Reproducability issue with svd_solver='arpack'
sc.pp.neighbors(S6, n_neighbors=10)
sc.tl.umap(S6)

sc.pl.pca_scatter(S6,color=['n_genes','n_counts','percent_mito','percent_ribo','phase','pt_group'],color_map=mymap)
sc.pl.umap(S6,color=['n_genes','n_counts','percent_mito','percent_ribo','phase','pt_group'],color_map=mymap)

sc.pp.pca(S8, n_comps=20, svd_solver='randomized') # Reproducability issue with svd_solver='arpack'
sc.pp.neighbors(S8, n_neighbors=10)
sc.tl.umap(S8)

sc.pl.pca_scatter(S8,color=['n_genes','n_counts','percent_mito','percent_ribo','phase','pt_group'],color_map=mymap)
sc.pl.umap(S8,color=['n_genes','n_counts','percent_mito','percent_ribo','phase','pt_group'],color_map=mymap)

## Estimate RNA velocity using dynamic model


scv.tl.recover_dynamics(S6)
scv.tl.velocity(S6, mode='dynamical')
scv.tl.velocity_graph(S6)

scv.tl.recover_dynamics(S8)
scv.tl.velocity(S8, mode='dynamical')
scv.tl.velocity_graph(S8)

# Plot the velocities in umap
scv.pl.velocity_embedding_stream(S6, basis='umap',color='pt_group',smooth=True,cutoff_perc=25,legend_loc='right marginal')
scv.pl.velocity_embedding_stream(S8, basis='umap',color='pt_group',smooth=True,cutoff_perc=25,legend_loc='right marginal')

# Paga graphs
S6.uns['neighbors']['distances'] = S6.obsp['distances']
S6.uns['neighbors']['connectivities'] = S6.obsp['connectivities']

scv.tl.paga(S6, groups='pt_group')
df = scv.get_df(S6, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

S8.uns['neighbors']['distances'] = S8.obsp['distances']
S8.uns['neighbors']['connectivities'] = S8.obsp['connectivities']

scv.tl.paga(S8, groups='pt_group')
df = scv.get_df(S8, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(S6, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5)
scv.pl.paga(S8, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5)

# noticed that arrow directions sometimes change in rerun, omitting arrows from final figs

#visualize gene sets
sc.pl.umap(S6, color=['S_score','G2M_score'],color_map=mymap)
sc.pl.umap(S8, color=['S_score','G2M_score'],color_map=mymap)
sc.pl.umap(S6, color=['S_score','G2M_score'],color_map=mymap)

#Make geneset list and name list using bash commands

gs_list = !ls gs_markers
name_list = [re.sub('_geneset.txt','',x) for x in gs_list]
gs_df = pd.DataFrame(gs_list)
name_df = pd.DataFrame(name_list)

#For-loop that reads in the lists and scores the geneset
for i in range(0,len(gs_df[0])):
    file_name = 'markers/'+gs_df[0][i]
    gs = pd.read_csv(file_name,header=None)
    gs = gs[np.isin(gs,adata.var_names)]
    gs = list(gs[0])
    score_name = name_df[0][i] + '_score'
    sc.tl.score_genes(adata,gene_list=gs,score_name=score_name)
    
adata_S6 = adata[adata.obs.batch == 'S6']
adata_S8 = adata[adata.obs.batch == 'S8']

for i in range(0,len(name_df[0])):
    obs_name = name_df[0][i] + '_score'
    S6.obs[obs_name] = adata_S6.obs[obs_name]
    
for i in range(0,len(name_df[0])):
    obs_name = name_df[0][i] + '_score'
    S8.obs[obs_name] = adata_S8.obs[obs_name]
    
sc.pl.umap(S6, color=[name_list[i] + '_score' for i in range(0,len(name_list))],color_map=mymap)
sc.pl.umap(S8, color=[name_list[i] + '_score' for i in range(0,len(name_list))],color_map=mymap)
