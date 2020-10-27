"""Rules related to viral progeny barcodes."""

rule viral_barcodes_in_progeny:
    """Parse and count viral barcodes from progeny sequencing data."""
    input:
        fastq_list=lambda wc: (expts.expt_viral_barcode_fastqs(wc.expt)
                               ['fastq_path']
                               .tolist()
                               ),
        viral_bc_locs=join(config['viral_fastq10x_dir'], 'viral_bc_locs.csv'),
        viral_genbank=config['viral_genbank'],
        notebook='notebooks/viral_barcodes_in_progeny.py.ipynb'
    output:
        viral_bc_in_progeny_csv=join(config['viral_progeny_dir'],
                                     "{expt}_viral_bc_in_progeny.csv.gz"),
    params:
        fastq_df=lambda wc: expts.expt_viral_barcode_fastqs(wc.expt),
        viral_barcode_length=config['viral_barcode_length'],
        viral_barcode_upstream_length=config['viral_barcode_upstream_length'],
        viral_barcode_mismatch_threshold=config['viral_barcode' \
                                                '_mismatch_threshold'],
        barcoded_viral_genes=barcoded_viral_genes
    log:
        notebook=join(config['log_dir'],
                      "{expt}_viral_barcodes_in_progeny.py.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/viral_barcodes_in_progeny.py.ipynb'
