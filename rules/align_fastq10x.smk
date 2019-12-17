"""``snakemake`` rules related aligned the 10X Illumina FASTQ reads."""


rule align_fastq10x_stats:
    """Aggregate 10X alignment stats and make summary plots."""
    input:
        summary=expand(join(config['aligned_fastq10x_dir'], "{sample10x}",
                            'Solo.out/Gene/Summary.csv'),
                       sample10x=samples_10x),
        umi_per_cell=expand(join(config['aligned_fastq10x_dir'], "{sample10x}",
                                 'Solo.out/Gene/UMIperCellSorted.txt'),
                            sample10x=samples_10x)
    output:
        stats=report(join(config['aligned_fastq10x_dir'], 'summary_stats.csv'),
                     caption='../report/align_fastq10x_stats.rst',
                     category='Aligning 10X FASTQs',
                     ),
        cells_plot=report(
                join(config['aligned_fastq10x_dir'], 'cells_plot.svg'),
                caption='../report/align_fastq10x_cells_plot.rst',
                category='Aligning 10X FASTQs',
                ),
        knee_plot=report(join(config['aligned_fastq10x_dir'], 'knee_plot.svg'),
                         caption='../report/align_fastq10x_knee_plot.rst',
                         category='Aligning 10X FASTQs',
                         ),
        per_cell_plot=report(
                join(config['aligned_fastq10x_dir'], 'per_cell_plot.svg'),
                caption='../report/align_fastq10x_per_cell_plot.rst',
                category='Aligning 10X FASTQs',
                ),
        mapping_rate_plot=report(
                join(config['aligned_fastq10x_dir'], 'mapping_rate_plot.svg'),
                caption='../report/align_fastq10x_mapping_rate_plot.rst',
                category='Aligning 10X FASTQs',
                )
    run:
        papermill.execute_notebook(
            input_path='notebooks/align_fastq10x_stats.ipynb',
            output_path=join(config['aligned_fastq10x_dir'],
                             'align_fastq10x_stats.ipynb'),
            cwd=os.getcwd(),
            parameters={
                'samples_10x': samples_10x,
                'input_summary': input.summary,
                'input_umi_per_cell': input.umi_per_cell,
                'output_stats': output.stats,
                'output_cells_plot': output.cells_plot,
                'output_knee_plot': output.knee_plot,
                'output_per_cell_plot': output.per_cell_plot,
                'output_mapping_rate_plot': output.mapping_rate_plot,
                },
            )


rule align_fastq10x:
    """Align 10X Illumina FASTQ reads with ``STARsolo``."""
    input:
        cb_whitelist_10x=config['cb_whitelist_10x'],
        refgenome=config['refgenome'],
        fastqR1=lambda wc: [join(config['fastq10x_dir'],
                                 f"{r}_all_R1.fastq.gz")
                            for r in (illumina_runs_10x
                                      .query(f"sample == '{wc.sample10x}'")
                                      .index
                                      )
                            ],
        fastqR2=lambda wc: [join(config['fastq10x_dir'],
                                 f"{r}_all_R2.fastq.gz")
                            for r in (illumina_runs_10x
                                      .query(f"sample == '{wc.sample10x}'")
                                      .index
                                      )
                            ]
    output:
        summary=join(config['aligned_fastq10x_dir'], "{sample10x}",
                     'Solo.out/Gene/Summary.csv'),
        umi_per_cell=join(config['aligned_fastq10x_dir'], "{sample10x}",
                          'Solo.out/Gene/UMIperCellSorted.txt'),
        matrix=join(config['aligned_fastq10x_dir'], "{sample10x}",
                    'Solo.out/Gene/filtered/matrix.mtx.gz'),
        features=join(config['aligned_fastq10x_dir'], "{sample10x}",
                      'Solo.out/Gene/filtered/features.tsv.gz'),
        barcodes=join(config['aligned_fastq10x_dir'], "{sample10x}",
                      'Solo.out/Gene/filtered/barcodes.tsv.gz'),
    params:
        outdir=join(config['aligned_fastq10x_dir'], "{sample10x}") + '/'
    threads: config['max_cpus']
    run:
        cmds = [
            'STAR',
            '--soloType', 'CB_UMI_Simple',
            '--genomeDir', input.refgenome,
            '--soloCBwhitelist', input.cb_whitelist_10x,
            '--readFilesIn', ','.join(input.fastqR2), ','.join(input.fastqR1),
            '--readFilesCommand', 'zcat',
            '--soloCBlen', str(config['cb_len_10x']),
            '--soloUMIlen', str(config['umi_len_10x']),
            '--soloUMIfiltering', 'MultiGeneUMI',
            '--soloCBmatchWLtype', '1MM_multi_pseudocounts',
            '--soloCellFilter', 'CellRanger2.2',
                                str(config['expect_ncells']), '0.99', '10',
            '--outSAMattributes', 'NH', 'HI', 'nM', 'AS', 'CB', 'UB',
            '--outSAMtype', 'BAM', 'SortedByCoordinate',
            '--runThreadN', str(threads),
            '--outFileNamePrefix', params.outdir,
            ]
        print(f"Running STARsolo with following command:\n{' '.join(cmds)}")
        os.makedirs(params.outdir, exist_ok=True)
        subprocess.check_call(cmds)

        # gzip filtered cell-gene matrix
        for fgz in [output.matrix, output.features, output.barcodes]:
            f = os.path.splitext(fgz)[0]
            print(f"gzipping {f}")
            assert os.path.isfile(f), f"Cannot find {f}"
            subprocess.check_call(['gzip', f])

rule get_cb_whitelist_10x:
    """Get whitelisted 10X cellbarcodes."""
    output: config['cb_whitelist_10x']
    params: ftp=config['cb_whitelist_10x_ftp']
    shell: "wget -O - {params.ftp} | gunzip -c > {output}"
