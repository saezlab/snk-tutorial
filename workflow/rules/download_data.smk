import os

checkpoint download:
    output:
        directory('data/filtered_gene_bc_matrices/hg19')
    shell:
        '(test -d data || mkdir data)'
        '&& wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz '
        '&& cd data; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz'

def get_downloaded_file(wildcards):
    checkpoint_output = checkpoints.download.get(**wildcards).output[0]
    return os.path.join(checkpoint_output, "matrix.mtx")

rule make_samples:
    input:
        get_downloaded_file
        # 'data/filtered_gene_bc_matrices/hg19/matrix.mtx'
    output:
        'data/sample1.h5ad',
        'data/sample2.h5ad',
        'data/sample3.h5ad'
    conda:
        '../envs/scanpy.yaml'
    script:
        '../scripts/fake_samples.py'
