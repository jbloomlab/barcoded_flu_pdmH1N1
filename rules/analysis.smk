"""``snakemake`` rules related to analysis of results."""


rule analyze_cell_gene_matrix:
    """Analyze the cell-gene matrix."""
    input:
        matrix=join(config['aligned_fastq10x_dir'], "{sample10x}",
                    'Solo.out/Gene/filtered/matrix.mtx'),
        features=join(config['aligned_fastq10x_dir'], "{sample10x}",
                      'Solo.out/Gene/filtered/features.tsv'),
        barcodes=join(config['aligned_fastq10x_dir'], "{sample10x}",
                      'Solo.out/Gene/filtered/barcodes.tsv'),
        viral_gtf=config['viral_gtf'],
        viraltag_counts=join(config['viral_fastq10x_dir'],
                             "viraltag_counts_{sample10x}.csv"),
        viralbc_counts=join(config['viral_fastq10x_dir'],
                            "viralbc_counts_{sample10x}.csv"),
        gapped_reads=join(config['viral_fastq10x_dir'], "gapped_reads.csv"),
        nb='notebooks/analyze_cell_gene_matrix.ipynb'
    output:
        nb=join(config['analysis_dir'],
                "{sample10x}_analyze_cell_gene_matrix.ipynb"),
        nb_html=report(join(config['analysis_dir'],
                            "{sample10x}_analyze_cell_gene_matrix.html"),
                       caption='../report/analyze_cell_gene_matrix.rst',
                       category='Analysis')
    run:
        run_nb_to_html(
                input_nb=input.nb,
                output_nb=output.nb,
                parameters={
                    'input_matrix': input.matrix,
                    'input_features': input.features,
                    'input_barcodes': input.barcodes,
                    'input_viral_gtf': input.viral_gtf,
                    'input_viraltag_counts': input.viraltag_counts,
                    'input_viralbc_counts': input.viralbc_counts,
                    'input_gapped_reads': input.gapped_reads,
                    },
                )
