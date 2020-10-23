"""Rules related to viral progeny barcodes."""

rule viral_barcodes_in_progeny:
    """Parse and count viral barcodes from progeny sequencing data."""
    input:
        fastq_list=lambda wc: expts.expt_viral_barcode_fastqs(wc.expt)['fastqs'].tolist(),
        viral_bc_locs=join(config['viral_fastq10x_dir'], 'viral_bc_locs.csv'),
        notebook='notebooks/viral_barcodes_in_progeny.py.ipynb'
    output:
        viral_bc_in_supernatant_csv=join(config['viral_progeny_dir'],
                                      "{expt}_viral_bc_in_supernatant.csv.gz"),
        viral_bc_in_secondinfection_csv=join(config['viral_progeny_dir'],
                                      "{expt}_viral_bc_in_secondinfection.csv.gz"),
    params:
        fastq_df=lambda wc: expts.expt_viral_barcode_fastqs(wc.expt),
        viral_barcode_upstream_length=config['viral_barcode_upstream_length'],
        barcoded_viral_genes=barcoded_viral_genes
    log:
        notebook=join(config['viral_progeny_dir'],
                      "{expt}_viral_barcodes_in_progeny.py.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/viral_barcodes_in_progeny.py.ipynb'