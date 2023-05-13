#Running RNA velocity analysis as presented in scvelo package documentation here for Nalm6 data
#Sample S2 = treatment, S3 = ctrl

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

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, dpi_save=300, vector_friendly=True)

# A nice color scheme for visualizing gene expression
colors_0 = plt.cm.OrRd(np.linspace(0.05, 1, 128))
colors_1 = plt.cm.Greys_r(np.linspace(0.8,0.9,20))
colors_Comb = np.vstack([colors_1, colors_0])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colors_Comb)

#Color palettes for plotting categorical variables 
#palette1 is Scanpy default and palette2 is for plotting pt_groups with consistency between samples

palette1 = ['#1f77b4',
            '#ff7f0e',
            '#279e68',
            '#d62728',
            '#aa40fc',
            '#8c564b',
            '#e377c2',
            '#b5bd61',
            '#17becf',
            '#aec7e8',
            '#ffbb78',
            '#98df8a',
            '#ff9896',
            '#c5b0d5',
            '#c49c94',
            '#f7b6d2',
            '#dbdb8d',
            '#9edae5',
            '#ad494a',
            '#8c6d31']

palette2 = ['#1f77b4',
            '#ff7f0e',
            '#279e68',
            '#d62728',
            #'#aa40fc',
            '#8c564b',
            '#e377c2',
            '#b5bd61',
            '#17becf',
            '#aec7e8',
            '#ffbb78',
            '#98df8a',
            '#ff9896',
            '#c5b0d5',
            '#c49c94',
            '#f7b6d2',
            '#dbdb8d',
            '#9edae5',
            '#ad494a',
            '#8c6d31']

# Read count data
S2 = read_10x_h5('cellranger/S2/outs/filtered_feature_bc_matrix.h5')
S3 = read_10x_h5('cellranger/S3/outs/filtered_feature_bc_matrix.h5')

S2s = scv.read('cellranger/S2/coding_velocyto_res/S2.loom')
S3s = scv.read('cellranger/S3/coding_velocyto_res/S3.loom')

# Get wanted barcodes and genes from text files to filter data
barcodes = pd.read_csv('cellbarcodes_pt_times_Nfilt.txt', sep="\t", header=None)
# Use the same variable genes as for RS4;11
vargenes = pd.read_csv('rs411_vargenes.txt',sep='\t',header=None)

# Change suffix to match the adata batches
barcodes[0] = [re.sub('-1_2','_S3',x) for x in barcodes[0]]
barcodes[0] = [re.sub('-1_1','_S2',x) for x in barcodes[0]]

pt_times = barcodes[1]
pt_times_list = list(pt_times)

# Edit cellIDs to match them
S2.obs_names = [re.sub('-1', '', x) for x in S2.obs_names]
S2s.obs_names = [re.sub('.*:|x', '', x) for x in S2s.obs_names]
S3.obs_names = [re.sub('-1', '', x) for x in S3.obs_names]
S3s.obs_names = [re.sub('.*:|x', '', x) for x in S3s.obs_names]

# Add spliced/unspliced counts as layers
S2 = scv.utils.merge(S2, S2s)
S3 = scv.utils.merge(S3, S3s)

# Combine datas
adata = S2.concatenate(S3, batch_categories=['S2', 'S3'], index_unique='_')

# Clean up AnnData object
scv.utils.cleanup(adata)
del adata.obs['Clusters'], adata.obs['_X'], adata.obs['_Y'], adata.obs['initial_size']

# Add pt-group to the observations
adata.obs['pt_group'] = pt_times_list

# Remove mouse genes
adata = adata[:,adata.var_names.str.startswith('hg19')]
adata.var_names = [re.sub('hg19_', '', x) for x in adata.var_names]

# Save raw counts
adata.raw = scv.pp.log1p(adata, copy=True)

scv.pl.proportions(adata)

# Prelim filter
sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=1)

mito_genes = [name for name in adata.var_names if name.startswith('MT-')]
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
ribo_genes = [name for name in adata.var_names if name.startswith('RP')]
adata.obs['percent_ribo'] = np.sum(adata[:, ribo_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

# Check QC metrics
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito', 'percent_ribo'], stripplot=False, multi_panel=True, groupby='batch')

scv.pp.normalize_per_cell(adata)
scv.pp.log1p(adata)

# Subdet adata to include only genes that are wanted to the analysis
adata = adata[:,np.isin(adata.var_names,vargenes[0])]

scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# Read cell cycle genes from gene list by Regev lab
s_genes = pd.read_csv('regev_lab_cell_cycle_genes_S_phase.txt', header=None)[0]
g2m_genes = pd.read_csv('regev_lab_cell_cycle_genes_G2M_phase.txt', header=None)[0]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

S2 = adata[adata.obs.batch == 'S2']

S3 = adata[adata.obs.batch == 'S3']

sc.pp.pca(S2, n_comps=20, svd_solver='randomized') # Reproducability issue with svd_solver='arpack'
sc.pp.neighbors(S2, n_neighbors=10)
sc.tl.umap(S2)

sc.pl.pca_scatter(S2,color=['n_genes','n_counts','percent_mito','percent_ribo','phase','pt_group'],color_map=mymap,palette=palette2,save='_pca_scatter_S2.png')
sc.pl.umap(S2,color=['n_genes','n_counts','percent_mito','percent_ribo','phase','pt_group'],color_map=mymap,palette=palette2,save='_umap_S2.png')

sc.pp.pca(S3, n_comps=20, svd_solver='randomized') # Reproducability issue with svd_solver='arpack'
sc.pp.neighbors(S3, n_neighbors=10)
sc.tl.umap(S3)

sc.pl.pca_scatter(S3,color=['n_genes','n_counts','percent_mito','percent_ribo','phase','pt_group'],color_map=mymap,save='_pca_scatter_S3.png')
sc.pl.umap(S3,color=['n_genes','n_counts','percent_mito','percent_ribo','phase','pt_group'],color_map=mymap,save='_umap_S3.png')

# Estimate RNA velocity, dynamic model

scv.tl.recover_dynamics(S2)
scv.tl.velocity(S2, mode='dynamical')
scv.tl.velocity_graph(S2)

scv.tl.recover_dynamics(S3)
scv.tl.velocity(S3, mode='dynamical')
scv.tl.velocity_graph(S3)

# Plot the velocities in umap
scv.pl.velocity_embedding_stream(S2, basis='umap',color='pt_group',legend_loc='right marginal',palette=palette2)
scv.pl.velocity_embedding_stream(S3, basis='umap',color='pt_group',legend_loc='right marginal')

# Plot the velocities in umap
scv.pl.velocity_embedding_stream(S2, basis='umap',color='pt_group',smooth=True, legend_loc='right marginal',palette=palette2,save='velocity_S2.png')
scv.pl.velocity_embedding_stream(S3, basis='umap',color='pt_group',smooth=True,legend_loc='right marginal',save='velocity_S3.png')

# Gene set scoring
sc.pl.umap(S2, color=['S_score','G2M_score'],color_map=mymap,save='_cellcyclescore_S2.png')
sc.pl.umap(S3, color=['S_score','G2M_score'],color_map=mymap,save='_cellcyclescore_S3.png')

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
    
adata_S2 = adata[adata.obs.batch == 'S2']
adata_S3 = adata[adata.obs.batch == 'S3']

for i in range(0,len(name_df[0])):
    obs_name = name_df[0][i] + '_score'
    S2.obs[obs_name] = adata_S2.obs[obs_name]
    
for i in range(0,len(name_df[0])):
    obs_name = name_df[0][i] + '_score'
    S3.obs[obs_name] = adata_S3.obs[obs_name]
    
sc.pl.umap(S2, color=[name_list[i] + '_score' for i in range(0,len(name_list))],color_map=mymap)
sc.pl.umap(S3, color=[name_list[i] + '_score' for i in range(0,len(name_list))],color_map=mymap)
