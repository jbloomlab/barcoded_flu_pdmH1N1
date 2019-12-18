"""``snakemake`` rules related analysis of results."""


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
        nb='notebooks/analyze_cell_gene_matrix.ipynb'
    output:
        nb=join(config['analysis_dir'],
                "{sample10x}_analyze_cell_gene_matrix.ipynb"),
        nb_html=report(join(config['analysis_dir'],
                            "{sample10x}_analyze_cell_gene_matrix.html"),
                       caption='../report/analyze_cell_gene_matrix.rst',
                       category='Analysis')
    run:
        papermill.execute_notebook(
            input_path=input.nb,
            output_path=output.nb,
            cwd=os.getcwd(),
            parameters={
                'input_matrix': input.matrix,
                'input_features': input.features,
                'input_barcodes': input.barcodes,
                'input_viral_gtf': input.viral_gtf,
                },
            )

        # https://github.com/ipython-contrib/jupyter_contrib_nbextensions/issues/901
        subprocess.check_call(['jupyter', 'nbconvert', output.nb,
                               '--to', 'html_embed', '--template', 'toc2'])
