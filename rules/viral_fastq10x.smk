"""Rules related viral reads in aligned 10x Illumina transcriptomics."""

rule viral_gene_progeny_relationship:
    """Plots a histogram of viral gene numbers binned by cells that contribute progeny."""
    input:
        viral_genes_by_cell_csv=join(config['viral_fastq10x_dir'],
                                     "{expt}_viral_genes_by_cell.csv.gz"),
        contributes_progeny_by_cell_csv=join(config['viral_fastq10x_dir'],
                                             ("{expt}_viral_gene_progeny"
                                              "_relationship.csv.gz")),
    output:
        plot=report(join(config['viral_fastq10x_dir'],
                         "{expt}_viral_gene_progeny_relationship.pdf"),
                    caption='../report/viral_gene_progeny_relationship.rst',
                    category="{expt}")
    log:
        notebook=join(config['log_dir'],
                      "viral_gene_progeny_relationship_{expt}.ipynb"
    conda: '../environment.yml'
    notebook:
        '../notebooks/viral_gene_progeny_relationship.py.ipynb'

rule contributes_progeny_by_cell:
    """Annotates each infected cell with whether it contributes progeny."""
    input:
        transcription_progeny_csv=join(config['viral_progeny_dir'],
                                       "{expt}_transcription_progeny.csv.gz"),
        viral_genes_by_cell_csv=join(config['viral_fastq10x_dir'],
                                     "{expt}_viral_genes_by_cell.csv.gz"),
    output:
        contributes_progeny_by_cell_csv=join(config['viral_fastq10x_dir'],
                                             ("{expt}_contributes_progeny"
                                              "_by_cell.csv.gz")),
        plot=report(join(config['viral_fastq10x_dir'],
                         "{expt}_contributes_progeny_by_cell.pdf"),
                    caption='../report/contributes_progeny_by_cell.rst',
                    category="{expt}")
    log:
        notebook=join(config['log_dir'],
                      "contributes_progeny_by_cell_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/contributes_progeny_by_cell.py.ipynb'

rule viral_gene_presence:
    """Call presence or absence of each viral gene in each cell."""
    input:
        matrix=join(config['aligned_fastq10x_dir'], "{expt}",
                    'Solo.out/GeneFull/filtered', 'matrix.mtx'),
        cell_barcodes=join(config['aligned_fastq10x_dir'], "{expt}",
                    'Solo.out/GeneFull/filtered', 'barcodes.tsv'),
        features=join(config['aligned_fastq10x_dir'], "{expt}",
                    'Solo.out/GeneFull/filtered', 'features.tsv'),
        cell_annotations=join(config['viral_tags_bcs_in_cells_dir'],
                              "{expt}_cell_barcodes_with_viral_tags.csv.gz"),
    output:
        viral_genes_by_cell_csv=join(config['viral_fastq10x_dir'],
                                     "{expt}_viral_genes_by_cell.csv.gz"),
        plot=report(join(config['viral_fastq10x_dir'],
                    "{expt}_viral_genes_by_cell.svg"),
                    caption='../report/viral_genes_by_cell.rst',
                    category="{expt}")
    params:
        viral_genes=viral_genes,
        barcoded_viral_genes=barcoded_viral_genes
    log:
        notebook=join(config['log_dir'],
                          "viral_gene_presence_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/viral_gene_presence.py.ipynb'

rule viral_transcript_coverage:
    """Coverage over viral transcripts in 10x transcriptomics."""
    input:
        bam=join(config['aligned_fastq10x_dir'], "{expt}",
                 'Aligned.sortedByCoord.out.bam'),
        bai=join(config['aligned_fastq10x_dir'], "{expt}",
                 'Aligned.sortedByCoord.out.bam.bai'),
        viral_genbank=config['viral_genbank'],
        notebook='notebooks/viral_transcript_coverage.py.ipynb'
    output:
        plot=report(join(config['viral_fastq10x_dir'],
                    "{expt}_viral_transcript_coverage.svg"),
                    caption='../report/viral_transcript_coverage.rst',
                    category="{expt}")
    log:
        notebook=join(config['log_dir'],
                      "viral_transcript_coverage_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/viral_transcript_coverage.py.ipynb'


rule viral_barcodes_by_cell_valid:
    """Calls valid viral barcodes in transcriptomic data."""
    input:
        viral_tag_by_cell_csv=join(config['viral_tags_bcs_in_cells_dir'],
                                   "{expt}_cell_barcodes_with_viral"
                                   "_tags.csv.gz"),
        viral_bc_by_cell_corrected_csv=join(config['viral_fastq10x_dir'],
                                            ("{expt}_viral_bc_by_cell"
                                             "_corrected.csv.gz")),
        viral_genes_by_cell_csv=join(config['viral_fastq10x_dir'],
                                     "{expt}_viral_genes_by_cell.csv.gz"),
    output:
        viral_barcodes_valid_csv=join(config['viral_fastq10x_dir'],
                                          ("{expt}_viral_bc_by"
                                           "_cell_valid.csv.gz")),
        plot=report(join(config['viral_fastq10x_dir'],
                         "{expt}_viral_barcodes_by_cell_valid.pdf"),
                    caption='../report/viral_bc_by_cell_valid.rst',
                    category="{expt}")
    params:
        barcoded_viral_genes=barcoded_viral_genes,
    log:
        notebook=join(config['log_dir'],
                      "viral_barcodes_by_cell_valid_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/viral_barcodes_by_cell_valid.py.ipynb'


rule correct_viral_barcodes_by_cell:
    """Correct viral barcodes within each cell."""
    input:
        viral_bc_by_cell_csv=join(config['viral_fastq10x_dir'],
                                  "{expt}_viral_bc_by_cell.csv.gz")
    output:
        viral_bc_by_cell_corrected_csv=join(config['viral_fastq10x_dir'],
                                            ("{expt}_viral_bc_by_cell_"
                                             "corrected.csv.gz")),
        plot=report(join(config['viral_fastq10x_dir'],
                         "{expt}_viral_bc_by_cell_corrected.pdf"),
                    caption='../report/viral_barcodes_by_cell_corrected.rst',
                    category="{expt}")
    log:
        notebook=join(config['log_dir'],
                      "correct_viral_barcodes_by_cell_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/correct_viral_barcodes_by_cell.py.ipynb'


rule viral_barcodes_by_cell:
    """Aggregate viral barcodes in transcriptomics by cell."""
    input:
        viral_bc_by_cell_umi_csv=join(config['viral_fastq10x_dir'],
                                      "{expt}_viral_bc_by_cell_umi.csv.gz"),
        notebook='notebooks/viral_barcodes_by_cell.py.ipynb'
    output:
        viral_bc_by_cell_csv=join(config['viral_fastq10x_dir'],
                                  "{expt}_viral_bc_by_cell.csv.gz"),
        plot=report(join(config['viral_fastq10x_dir'],
                         "{expt}_viral_bc_by_cell.svg"),
                    caption='../report/viral_barcodes_by_cell.rst',
                    category="{expt}")
    log:
        notebook=join(config['log_dir'],
                      "viral_barcodes_by_cell_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/viral_barcodes_by_cell.py.ipynb'


rule viral_tags_by_cell:
    """Aggregate viral tags in transcriptomics by cell."""
    input:
        viral_tag_by_cell_umi_csv=join(config['viral_fastq10x_dir'],
                                       "{expt}_viral_tag_by_cell_umi.csv.gz"),
        notebook='notebooks/viral_tags_by_cell.py.ipynb'
    output:
        viral_tag_by_cell_csv=join(config['viral_fastq10x_dir'],
                                   "{expt}_viral_tag_by_cell.csv.gz"),
        plot=report(join(config['viral_fastq10x_dir'],
                         "{expt}_viral_tag_by_cell.svg"),
                    caption='../report/viral_tags_by_cell.rst',
                    category="{expt}")
    log:
        notebook=join(config['log_dir'],
                      "viral_tags_by_cell_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/viral_tags_by_cell.py.ipynb'


rule viral_barcodes_in_transcripts:
    """Extract viral barcodes from 10x transcriptomic alignments."""
    input:
        bam=join(config['aligned_fastq10x_dir'], "{expt}",
                 'Aligned.sortedByCoord.out.bam'),
        bai=join(config['aligned_fastq10x_dir'], "{expt}",
                 'Aligned.sortedByCoord.out.bam.bai'),
        cell_barcodes=join(config['aligned_fastq10x_dir'], "{expt}",
                           'Solo.out/GeneFull/filtered/barcodes.tsv'),
        viral_bc_locs=join(config['viral_fastq10x_dir'], 'viral_bc_locs.csv'),
        notebook='notebooks/viral_barcodes_in_transcripts.py.ipynb'
    output:
        viral_bc_by_cell_umi_csv=join(config['viral_fastq10x_dir'],
                                      "{expt}_viral_bc_by_cell_umi.csv.gz"),
    log:
        notebook=join(config['log_dir'],
                      "viral_barcodes_in_transcripts_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/viral_barcodes_in_transcripts.py.ipynb'


rule viral_tags_in_transcripts:
    """Extract viral tags from 10x transcriptomic alignments."""
    input:
        bam=join(config['aligned_fastq10x_dir'], "{expt}",
                 'Aligned.sortedByCoord.out.bam'),
        bai=join(config['aligned_fastq10x_dir'], "{expt}",
                 'Aligned.sortedByCoord.out.bam.bai'),
        cell_barcodes=join(config['aligned_fastq10x_dir'], "{expt}",
                           'Solo.out/GeneFull/filtered/barcodes.tsv'),
        viral_tag_locs=join(config['viral_fastq10x_dir'],
                            'viral_tag_locs.csv'),
        viral_tag_identities=config['viral_tag_identities'],
        notebook='notebooks/viral_tags_in_transcripts.py.ipynb'
    output:
        viral_tag_by_cell_umi_csv=join(config['viral_fastq10x_dir'],
                                       "{expt}_viral_tag_by_cell_umi.csv.gz"),
    log:
        notebook=join(config['log_dir'],
                      "viral_tags_in_transcripts_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/viral_tags_in_transcripts.py.ipynb'


rule viral_bc_locs:
    """Locations of viral barcodes in 1-based indexing."""
    input:
        viral_genbank=config['viral_genbank']
    output:
        viral_bc_locs=join(config['viral_fastq10x_dir'], 'viral_bc_locs.csv'),
    conda: '../environment.yml'
    log:
        log=join(config['log_dir'], 'viral_bc_locs.log')
    script:
        '../scripts/viral_bc_locs.py'


rule viral_tag_locs:
    """Locations of viral tags in 1-based indexing."""
    input:
        viral_genbank=config['viral_genbank']
    output:
        viral_tag_locs=join(config['viral_fastq10x_dir'],
                            'viral_tag_locs.csv'),
    conda: '../environment.yml'
    log:
        log=join(config['log_dir'], 'viral_tag_locs.log')
    script:
        '../scripts/viral_tag_locs.py'
