#!/usr/bin/env python3
import numpy as np, pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import phate
import umap

#%%
YEOH_RPATH = "data/GSE67684/processed/subset_yeoh.tsv"
data = pd.read_csv(YEOH_RPATH, sep="\t")

META_RPATH = "data/GSE67684/processed/full_metadata.tsv"
full_metadata = pd.read_csv(META_RPATH, sep="\t", dtype="category")
full_metadata.batch_info = full_metadata.batch_info.astype("int16")

MNN_YEOH_RPATH = "dump/yeoh_mnn.tsv"
mnn_data = pd.read_csv(MNN_YEOH_RPATH, sep="\t")

#%%
### METADATA
metadata = full_metadata.loc[list(data)]
# metadata["col"] = metadata.class_info.str.cat(metadata.label).astype("category")
# X = data_pca
# fig1, ax1 = plt.subplots(1, 1, figsize=(12,4))
# # Matplotlib only plots integers as col
# # ax1.scatter(X[:,0], X[:,1],
# #             c=metadata.batch_info)
# # Seaborn plots categorical variables as col
# sns.scatterplot(X[:,0], X[:,1],
#                 hue=batch_info,
#                 linewidth=0, s=60, legend=False, ax=ax1)

# pal = sns.color_palette(["skyblue","blue","orange","red","green"])
#%%
def plotEval(X, metadata, main):
    # Error when category dtype is able to be coerced to float
    # Solution: Add alphabet in front of batch number
    # batch_info = metadata.batch_info.apply(lambda x: "B"+x)
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15,10))
    sns.scatterplot(X[:,0], X[:,1],
                    hue=metadata.batch_info, palette="viridis",
                    linewidth=0, s=60, legend=False, ax=ax1)
    sns.scatterplot(X[:,0], X[:,1],
                    hue=metadata.class_info,
                    linewidth=0, s=60, legend=False, ax=ax2)
    sns.scatterplot(X[:,0], X[:,1],
                    hue=metadata.subtype,
                    linewidth=0, s=60, legend=False, ax=ax3)
    sns.scatterplot(X[:,0], X[:,1],
                    hue=metadata.label,
                    linewidth=0, s=60, legend=False, ax=ax4)
    ax1.set_title("Batch")
    ax2.set_title("Timepoint")
    ax3.set_title("Subtype")
    ax4.set_title("Relapse")
    fig.savefig("dump/original_{}.pdf".format(main), bbox_inches="tight")

#%%
pca = PCA(n_components = 2)
data_pca = pca.fit_transform(data.T) # Remember to transpose data!

plotEval(data_pca, metadata, "pca")

#%%
tsne2D = TSNE(n_components=2, perplexity=30, random_state=1)
data_tsne2D = tsne2D.fit_transform(data.T)

plotEval(data_tsne2D, metadata, "tsne")
#%%
phate_op = phate.PHATE()
# PHATE: 3D
# phate_operator.set_params(n_components=3)
data_phate = phate_op.fit_transform(data.T)

plotEval(data_phate, metadata, "phate")

#%%
umap_op = umap.UMAP(n_neighbors=15)
data_umap = umap_op.fit_transform(data.T)

plotEval(data_umap, metadata, "umap")


# Check batch timepoints. Change colour?
# Plot dates as colours instead?

# Try different parameters for the DR methods?
