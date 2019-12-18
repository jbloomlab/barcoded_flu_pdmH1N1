"""``snakemake`` rules related aligned the 10X Illumina FASTQ reads."""


rule align_fastq10x_summary:
    """Summarize 10X FASTQ alignments."""
    input:
        summary=expand(join(config['aligned_fastq10x_dir'], "{sample10x}",
                            'Solo.out/Gene/Summary.csv'),
                       sample10x=samples_10x),
        umi_per_cell=expand(join(config['aligned_fastq10x_dir'], "{sample10x}",
                                 'Solo.out/Gene/UMIperCellSorted.txt'),
                            sample10x=samples_10x),
        nb='notebooks/align_fastq10x_summary.ipynb'
    output:
        nb=join(config['aligned_fastq10x_dir'],
                'align_fastq10x_summary.ipynb'),
        nb_html=report(join(config['aligned_fastq10x_dir'],
                            'align_fastq10x_summary.html'),
                       caption='../report/align_fastq10x_summary.rst',
                       category='Aligning 10X FASTQs')
    run:
        papermill.execute_notebook(
            input_path=input.nb,
            output_path=output.nb,
            cwd=os.getcwd(),
            parameters={
                'samples_10x': samples_10x,
                'input_summary': input.summary,
                'input_umi_per_cell': input.umi_per_cell,
                },
            )

        # https://github.com/ipython-contrib/jupyter_contrib_nbextensions/issues/901
        subprocess.check_call(['jupyter', 'nbconvert', output.nb,
                               '--to', 'html_embed', '--template', 'toc2'])


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
                    'Solo.out/Gene/filtered/matrix.mtx'),
        features=join(config['aligned_fastq10x_dir'], "{sample10x}",
                      'Solo.out/Gene/filtered/features.tsv'),
        barcodes=join(config['aligned_fastq10x_dir'], "{sample10x}",
                      'Solo.out/Gene/filtered/barcodes.tsv'),
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


rule get_cb_whitelist_10x:
    """Get whitelisted 10X cellbarcodes."""
    output: config['cb_whitelist_10x']
    params: ftp=config['cb_whitelist_10x_ftp']
    shell: "wget -O - {params.ftp} | gunzip -c > {output}"
