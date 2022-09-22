# %%
import scanpy as sc
import logging

print('scanpy version {0}'.format(sc.__version__))

import numpy as np
import pandas as pd

# %%

if 'snakemake' in locals():
    filepath = snakemake.input[0]

    if len(snakemake.params) == 0:
        QC_params = {'min_gene': 200, 'min_cells': 3, 'max_pct_mt': 5}
        logging.warning('No QC parameters were provided! Using defaults of min_gene = {0}, min_cells = {1} and max_pct_mt = {2}'.format(QC_params['min_gene'], QC_params['min_cells'], QC_params['max_pct_mt']))
    
    elif len(snakemake.params) == 3:
        QC_params = {'min_gene': snakemake.params[0], 'min_cells': snakemake.params[1], 'max_pct_mt': snakemake.params[2]}
    else:
        QC_params = snakemake.params[0]
    
else:
    filepath = '../../data/sample1.h5ad'
    QC_params = {'min_gene': 200, 'min_cells': 3, 'max_pct_mt': 5}

# %%
adata = sc.read_h5ad(filepath)

print('The raw data contains {0} cells and {1} genes'.format(adata.n_obs, adata.n_vars))

# %%
# Basic filtering
sc.pp.filter_cells(adata, min_genes = QC_params['min_gene'])
sc.pp.filter_genes(adata, min_cells = QC_params['min_cells'])

# Compute QC metrics (mice mitochondrial genes start with mt-, not MT- (humans))
adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < QC_params['max_pct_mt'], :]

# %%
sc.pp.normalize_total(adata, target_sum=1e4)

# %%
if 'snakemake' in locals():
    adata.write_h5ad(snakemake.output[0])
    print('The filtered data contains {0} cells and {1} genes'.format(adata.n_obs, adata.n_vars))
    print(adata)
