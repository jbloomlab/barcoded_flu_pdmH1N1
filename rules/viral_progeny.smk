"""Rules related to viral progeny barcodes."""

rule process_viral_barcode_replicates:
    """Plot correlation and average viral barcode replicates."""
    input:
        viral_bc_in_progeny_csv=join(config['viral_progeny_dir'],
                                     "{expt}_viral_bc_in_progeny.csv.gz"),
    output:
        viral_bc_in_progeny_freq_csv=join(config['viral_progeny_dir'],
                                          "{expt}_"
                                          "viral_bc_in_progeny_freq.csv.gz"),

rule viral_barcodes_in_progeny:
    """Parse and count viral barcodes from progeny sequencing data."""
    input:
        fastq_list=lambda wc: (expts.expt_viral_barcode_fastqs(wc.expt)
                               ['fastq_path']
                               .tolist()
                               ),
        viral_genbank=config['viral_genbank'],
        notebook='notebooks/viral_barcodes_in_progeny.py.ipynb'
    output:
        viral_bc_in_progeny_csv=join(config['viral_progeny_dir'],
                                     "{expt}_viral_bc_in_progeny.csv.gz"),
        viral_bc_fates_csv=join(config['viral_progeny_dir'],
                                "{expt}_viral_bc_fates.csv.gz"),
        plot=report(join(config['viral_progeny_dir'],
                         "{expt}_viral_bc_fates.svg"),
                    caption='../report/viral_bc_fates.rst',
                    category="{expt}")
    params:
        fastq_df=lambda wc: expts.expt_viral_barcode_fastqs(wc.expt),
        viral_barcode_upstream_length=config['viral_barcode_upstream_length'],
        viral_barcode_mismatch=config['viral_barcode_mismatch'],
        barcoded_viral_genes=barcoded_viral_genes
    log:
        notebook=join(config['log_dir'],
                      "viral_barcodes_in_progeny_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/viral_barcodes_in_progeny.py.ipynb'
