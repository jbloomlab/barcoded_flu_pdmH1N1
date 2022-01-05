"""Rules related to viral progeny barcodes."""

rule transcription_progeny_correlation:
    """Plot the correlation between viral transcription and progeny production."""
    input:
        cell_annotations=join(config['viral_tags_bcs_in_cells_dir'],
                              "{expt}_cell_barcodes_with_viral_tags.csv.gz"),
        viral_barcodes_valid_csv=join(config['viral_fastq10x_dir'],
                                          ("{expt}_viral"
                                           "_bc_by_cell_valid.csv.gz")),
        filtered_progeny_viral_bc_csv=join(config['viral_progeny_dir'],
                                          "{expt}_"
                                          "filtered_progeny_viral_bc.csv.gz"),
    output:
        transcription_progeny_csv=join(config['viral_progeny_dir'],
                                       "{expt}_transcription_progeny.csv.gz"),
        plot=report(join(config['viral_progeny_dir'],
                         "{expt}_transcription_progeny_correlation.pdf"),
                         caption='../report/transcription_progeny_correlation.rst',
                         category="{expt}")
    params:
        viral_genes=viral_genes,
        barcoded_viral_genes=barcoded_viral_genes
    log:
        notebook=join(config['log_dir'],
                      "transcription_progeny_correlation_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/transcription_progeny_correlation.py.ipynb'

rule filter_progeny:
    """Filter and average viral barcode replicates in progeny."""
    input:
        cell_annotations=join(config['viral_tags_bcs_in_cells_dir'],
                              "{expt}_cell_barcodes_with_viral_tags.csv.gz"),
        viral_bc_in_progeny_corrected_csv=join(config['viral_progeny_dir'],
                                               ("{expt}_viral_bc_in_progeny_"
                                                "corrected.csv.gz")),
        viral_barcodes_valid_csv=join(config['viral_fastq10x_dir'],
                                          ("{expt}_viral"
                                           "_bc_by_cell_valid.csv.gz")),
    output:
        filtered_progeny_viral_bc_csv=join(config['viral_progeny_dir'],
                                          "{expt}_"
                                          "filtered_progeny_viral_bc.csv.gz"),
        plot=report(join(config['viral_progeny_dir'],
                         "{expt}_filtered_progeny_viral_bc.pdf"),
                    caption='../report/filter_progeny.rst',
                    category="{expt}")
    log:
        notebook=join(config['log_dir'],
                      "filter_progeny_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/filter_progeny.py.ipynb'

rule correct_viral_barcodes_in_progeny:
    """Correct viral barcodes in progeny."""
    input:
        viral_bc_in_progeny_csv=join(config['viral_progeny_dir'],
                                     "{expt}_viral_bc_in_progeny.csv.gz"),
    output:
        viral_bc_in_progeny_corrected_csv=join(config['viral_progeny_dir'],
                                               ("{expt}_viral_bc_in_progeny_"
                                                "corrected.csv.gz")),
        plot=report(join(config['viral_progeny_dir'],
                         "{expt}_viral_bc_in_progeny_corrected.pdf"),
                    caption='../report/viral_barcodes_in_progeny_corrected.rst',
                    category="{expt}")
    log:
        notebook=join(config['log_dir'],
                      "correct_viral_barcodes_in_progeny_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/correct_viral_barcodes_in_progeny.py.ipynb'

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
        viral_barcode_minq=config['viral_barcode_minq'],
        barcoded_viral_genes=barcoded_viral_genes
    log:
        notebook=join(config['log_dir'],
                      "viral_barcodes_in_progeny_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/viral_barcodes_in_progeny.py.ipynb'
