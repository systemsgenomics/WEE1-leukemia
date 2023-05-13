#RS4;11 10 h scRNA-seq velocity analysis with dynamical model
#Running RNA velocity analysis as presented in scvelo package documentation here for RS4;11 data
#Sample M2 = treatment, M1 = ctrl

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
            '#c7c7c7',
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

adata = sc.read('RS411_10h.h5ad')
adata.obs['pt_group_cutoff0.5'] = adata.obs['pt_group']
df = pd.concat([meta_M1,meta_M2],axis=0)
adata.obs['pt_group'] = df['predicted.id']

# Read count data
M1 = read_10x_h5('M1R-SCI7T014-SCI5T014_HK2L7DSX3/outs/filtered_feature_bc_matrix.h5')
M2 = read_10x_h5('M2R-SCI7T026-SCI5T026_HK2L7DSX3/outs/filtered_feature_bc_matrix.h5')

# Read spliced/unspliced count data
M1s = scv.read('velocyto/M1R-SCI7T014-SCI5T014_HK2L7DSX3.loom')
M2s = scv.read('velocyto/split1_MOAU5.loom')

# Get wanted barcodes from text file to filter data
barcodes_M1 = pd.read_csv('cellbarcodes_RS411_DMSO.10h.txt', sep="\t", header=0)
barcodes_M2 = pd.read_csv('cellbarcodes_RS411_Wee1i.10h.txt', sep="\t", header=0)

meta_M1 = pd.read_csv('meta.RS41124h.geneActivity.withPredictedRNAcluLabels_RS411_10h.DMSO.txt', sep="\t", header=0)
meta_M2 = pd.read_csv('meta.RS41124h.geneActivity.withPredictedRNAcluLabels_RS411_10h.wee1i.txt', sep="\t", header=0)

# Change suffix to match the adata batches
barcodes_M1['x'] = [re.sub('-1_1','',x) for x in barcodes_M1['x']]
barcodes_M2['x'] = [re.sub('-1_2','',x) for x in barcodes_M2['x']]

# Change suffix to match the adata batches
meta_M1.index = [re.sub('-1_1','_M1',x) for x in meta_M1.index]
meta_M2.index = [re.sub('-1_2','_M2',x) for x in meta_M2.index]

# Add spliced/unspliced counts as layers
M1 = scv.utils.merge(M1, M1s)
M2 = scv.utils.merge(M2,M2s)

M1 = M1[np.isin(M1.obs_names,barcodes_M1['x'])]
M2 = M2[np.isin(M2.obs_names,barcodes_M2['x'])]

M1.obs['pt_group'] = meta_M1['predicted.idv2']
M2.obs['pt_group'] = meta_M2['predicted.idv2']

# Combine datas
adata = M1.concatenate(M2, batch_categories=['M1', 'M2'], index_unique='_')

# Clean up AnnData object
scv.utils.cleanup(adata)

# Save raw counts
adata.raw = scv.pp.log1p(adata, copy=True)

scv.pl.proportions(adata,groupby='batch')

# Prelim filter
sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=10)

mito_genes = [name for name in adata.var_names if name.startswith('MT-')]
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
ribo_genes = [name for name in adata.var_names if name.startswith('RP')]
adata.obs['percent_ribo'] = np.sum(adata[:, ribo_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

# Check QC metrics
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito', 'percent_ribo'], stripplot=False, multi_panel=True,groupby='batch')

scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)

scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# Read cell cycle genes from gene list by Regev lab
s_genes = pd.read_csv('regev_lab_cell_cycle_genes_S_phase.txt', header=None)[0]
g2m_genes = pd.read_csv('regev_lab_cell_cycle_genes_G2M_phase.txt', header=None)[0]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)


M1 = adata[adata.obs.batch == 'M1']
M2 = adata[adata.obs.batch == 'M2']

sc.tl.leiden(M1,resolution=1,key_added='leiden_r1')
sc.tl.leiden(M2,resolution=1,key_added='leiden_r1')

sc.pp.pca(M1, n_comps=20, svd_solver='randomized') # Reproducability issue with svd_solver='arpack'
sc.pp.neighbors(M1, n_neighbors=10)
sc.tl.umap(M1)

sc.pl.pca_scatter(M1,color=['n_genes','n_counts','percent_mito','percent_ribo','phase','pt_group'],color_map=mymap,palette=palette1)
sc.pl.umap(M1,color=['n_genes','n_counts','percent_mito','percent_ribo','phase','pt_group'],color_map=mymap,palette=palette1)

sc.pp.pca(M2, n_comps=20, svd_solver='randomized') # Reproducability issue with svd_solver='arpack'
sc.pp.neighbors(M2, n_neighbors=10)
sc.tl.umap(M2)

sc.pl.pca_scatter(M2,color=['n_genes','n_counts','percent_mito','percent_ribo','phase','pt_group'],color_map=mymap,palette=palette1)
sc.pl.umap(M2,color=['n_genes','n_counts','percent_mito','percent_ribo','phase','pt_group'],color_map=mymap,palette=palette1)

# Estimate RNA velocities, dynamic model
scv.tl.recover_dynamics(M1)
scv.tl.velocity(M1, mode='dynamical')
scv.tl.velocity_graph(M1)

scv.tl.recover_dynamics(M2)
scv.tl.velocity(M2, mode='dynamical')
scv.tl.velocity_graph(M2)

# Paga graphs
M1.uns['neighbors']['distances'] = M1.obsp['distances']
M1.uns['neighbors']['connectivities'] = M1.obsp['connectivities']

scv.tl.paga(M1, groups='pt_group')
df = scv.get_df(M1, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

M2.uns['neighbors']['distances'] = M2.obsp['distances']
M2.uns['neighbors']['connectivities'] = M2.obsp['connectivities']

scv.tl.paga(M2, groups='pt_group')
df = scv.get_df(M2, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(M1, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5,palette=palette1)
scv.pl.paga(M2, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5,palette=palette1)

sc.pl.umap(M1, color=['S_score','G2M_score'],color_map=mymap)
sc.pl.umap(M2, color=['S_score','G2M_score'],color_map=mymap)

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
    
adata_M1 = adata[adata.obs.batch == 'M1']
adata_M2 = adata[adata.obs.batch == 'M2']

for i in range(0,len(name_df[0])):
    obs_name = name_df[0][i] + '_score'
    M1.obs[obs_name] = adata_M1.obs[obs_name]
    
for i in range(0,len(name_df[0])):
    obs_name = name_df[0][i] + '_score'
    M2.obs[obs_name] = adata_M2.obs[obs_name]
    
sc.pl.umap(M1, color=[name_list[i] + '_score' for i in range(0,len(name_list))],color_map=mymap)
sc.pl.umap(M2, color=[name_list[i] + '_score' for i in range(0,len(name_list))],color_map=mymap)
