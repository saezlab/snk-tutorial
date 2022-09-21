# %%
import scanpy as sc
import os

# %%
if 'snakemake' in locals():
    dp = os.path.dirname(snakemake.input[0])
else:
    dp = '../../data/filtered_gene_bc_matrices/hg19/'

# %%
adata = sc.read_10x_mtx(path = dp)
adata.var_names_make_unique()
print(adata)

# %%
if 'snakemake' in locals():
    adata[0:900].write(snakemake.output[0])
    adata[900:1800].write(snakemake.output[1])
    adata[1800:2700].write(snakemake.output[2])

