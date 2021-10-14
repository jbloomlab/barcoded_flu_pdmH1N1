"""Rules annotating cells by viral tags / barcodes in transcriptomic data."""


rule assign_viral_tags_to_cells:
    """Assign infected cells their viral tag variants."""
    input:
        viral_tag_by_cell_csv=join(config['viral_fastq10x_dir'],
                                   "{expt}_viral_tag_by_cell.csv.gz"),
        infection_status_csv=join(config['viral_fastq10x_dir'],
                                  "{expt}_infection_status_csv.gz"),
        notebook='notebooks/assign_viral_tags_to_cells.py.ipynb'
    output:
        cell_annotations=join(config['viral_tags_bcs_in_cells_dir'],
                              "{expt}_cell_barcodes_with_viral_tags.csv.gz"),
        plot=report(join(config['viral_tags_bcs_in_cells_dir'],
                         "{expt}_assign_viral_tags_to_cells.svg"),
                    caption='../report/assign_viral_tags_to_cells.rst',
                    category="{expt}")
    params:
        viral_genes=viral_genes,
        viral_tags=viral_tags,
        fdr=config['viral_tag_fdr'],
    log:
        notebook=join(config['log_dir'],
                      "assign_viral_tags_to_cells_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/assign_viral_tags_to_cells.py.ipynb'

rule assign_infection_status:
    """Identify infected cells by viral burden."""
    input:
        matrix=join(config['aligned_fastq10x_dir'], "{expt}",
                    'Solo.out/GeneFull/filtered/matrix.mtx'),
        features=join(config['aligned_fastq10x_dir'], "{expt}",
                      'Solo.out/GeneFull/filtered/features.tsv'),
        cell_barcodes=join(config['aligned_fastq10x_dir'], "{expt}",
                           'Solo.out/GeneFull/filtered/barcodes.tsv'),
        notebook='notebooks/assign_infection_status.py.ipynb'
    output:
        infection_status_csv=join(config['viral_fastq10x_dir'],
                                  "{expt}_infection_status_csv.gz"),
        plot=report(join(config['viral_tags_bcs_in_cells_dir'],
                         "{expt}_assign_infection_status.svg"),
                    caption='../report/assign_infection_status.rst',
                    category="{expt}")
    params:
        viral_genes=viral_genes,
        infection_threshold=config['infection_threshold'],
    log:
        notebook=join(config['log_dir'],
                      "assign_infection_status_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/assign_infection_status.py.ipynb'
