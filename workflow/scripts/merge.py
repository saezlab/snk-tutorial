# %%
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# %%
if 'snakemake' in locals():
    filepaths = snakemake.input

else:
    filepaths = ['../../results/QC/filtered_sample1.h5ad', '../../results/QC/filtered_sample2.h5ad']


# %%
print('Loading datasets')
# combine sample data
adatas = []
for i in range(len(filepaths)):
    print('Loading {0}'.format(filepaths[i][0:-5].split('_')[1]))
    adata = sc.read_h5ad(filepaths[i])
    adatas.append(adata)

# Merge objects and delete list
adata = adatas[0].concatenate(adatas[1:], join='outer')
del adatas

#clean adata (-> else cannot write to file)
mt_columns = [x.startswith('mt') for x in adata.var.columns.to_list()]
col = adata.var.columns[mt_columns].to_list()
adata.var = adata.var.drop(col, axis = 1)
adata.var['mt'] = adata.var_names.str.startswith('mt-')

print('Merged samples')
print(adata)

# %%
print('Log-transforming dataset')
sc.pp.log1p(adata)

# %%
print('Selecting highly variable genes')
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

# %%
print('Computing principle components')
sc.tl.pca(adata, svd_solver='arpack')

# %%
print('Computing embeddings')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
print('Leiden clustering')
sc.tl.leiden(adata)

# %%
if 'snakemake' in locals():

    #clean adata.var
    #remove duplicate columns
    gene_index = pd.Series(adata.var.columns).str.contains('gene_ids').tolist()
    genes = adata.var.loc[:,gene_index]
    temp = adata.var.loc[:,np.invert(gene_index)]

    type_index = pd.Series(temp.columns).str.contains('feature_types').tolist()
    exp_type = temp.loc[:,type_index]
    temp = temp.loc[:,np.invert(type_index)]

    temp["gene_ids"] = [row[pd.Series(row).str.contains('ENS').tolist().index(True)] for _, row in genes.iterrows()]


    #fill nan with 0s
    adata.var = temp.fillna(0)
    adata.var = adata.var.drop(["mt","highly_variable"],axis=1)

    adata.write(snakemake.output[0])
    
    print('Saved data to file: {0}'.format(snakemake.output[0]))

# %%
fig, axes = plt.subplots(1,3,figsize = (15, 5))
axes = axes.flatten()
fig.suptitle('UMAP embedding', fontsize=14)

sc.pl.umap(adata, color='leiden', ax = axes[0], show = False)
sc.pl.umap(adata, color='CST3', ax = axes[1], show = False)
sc.pl.umap(adata, color='NKG7', ax = axes[2], show = False)

fig.set_facecolor('white')
plt.tight_layout()

if 'snakemake' in locals():
    print('Saved plot to file: {0}'.format(snakemake.output[1]))
    plt.savefig(snakemake.output[1], dpi=300)


# %%



