"""Rules to annotate cells by viral tags / barcodes from transcriptomic data."""


rule assign_viral_tags_to_cells:
    """Identify infected cells and assign them their viral tag variants."""
    input:
        viral_tag_by_cell_csv=join(config['viral_fastq10x_dir'],
                                   "{expt}_viral_tag_by_cell.csv.gz"),
        matrix=join(config['aligned_fastq10x_dir'], "{expt}",
                    'Solo.out/GeneFull/filtered/matrix.mtx'),
        features=join(config['aligned_fastq10x_dir'], "{expt}",
                      'Solo.out/GeneFull/filtered/features.tsv'),
        cell_barcodes=join(config['aligned_fastq10x_dir'], "{expt}",
                           'Solo.out/GeneFull/filtered/barcodes.tsv'),
        notebook='notebooks/assign_viral_tags_to_cells.py.ipynb'
    output:
        cell_annotations=join(config['viral_tags_bcs_in_cells_dir'],
                              "{expt}_cell_barcodes_with_viral_tags.csv"),
        plot=report(join(config['viral_tags_bcs_in_cells_dir'],
                         "{expt}_assign_viral_tags_to_cells.svg"),
                    caption='../report/assign_viral_tags_to_cells.rst',
                    category="{expt}")
    params:
        viral_genes=viral_genes,
        viral_tags=viral_tags,
        fdr=config['infection_calling_by_viral_tag_fdr'],
    log:
        notebook=join(config['viral_tags_bcs_in_cells_dir'],
                      "{expt}_assign_viral_tags_to_cells.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/assign_viral_tags_to_cells.py.ipynb'
