"""``snakemake`` rules related aligned the 10X Illumina FASTQ reads."""


rule fastq10x_transcript_coverage:
    """Summarize transcript coverage for a few selected genes."""
    input:
        fastq10x_bams=expand(join(config['aligned_fastq10x_dir'],
                                  "{sample10x}/Aligned.sortedByCoord.out.bam"),
                              sample10x=samples_10x),
        fastq10x_bais=expand(join(config['aligned_fastq10x_dir'],
                                  "{sample10x}",
                                  'Aligned.sortedByCoord.out.bam.bai'),
                             sample10x=samples_10x),
        viral_gtf=config['viral_gtf'],
        gtf=join(config['genome_dir'], 'cell_and_virus_gtf.gtf'),
        nb='notebooks/fastq10x_transcript_coverage.ipynb'
    output:
        nb=join(config['aligned_fastq10x_dir'],
                'fastq10x_transcript_coverage.ipynb'),
        nb_html=report(join(config['aligned_fastq10x_dir'],
                            'fastq10x_transcript_coverage.html'),
                       caption='../report/fastq10x_transcript_coverage.rst',
                       category='Aligning 10X FASTQs')
    run:
        run_nb_to_html(input_nb=input.nb,
                       output_nb=output.nb,
                       parameters={
                            'samples_10x': samples_10x,
                            'input_fastq10x_bams': input.fastq10x_bams,
                            'input_fastq10x_bais': input.fastq10x_bais,
                            'input_viral_gtf': input.viral_gtf,
                            'input_gtf': input.gtf,
                            },
                       )


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
        run_nb_to_html(input_nb=input.nb,
                       output_nb=output.nb,
                       parameters={
                            'samples_10x': samples_10x,
                            'input_summary': input.summary,
                            'input_umi_per_cell': input.umi_per_cell,
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
                    'Solo.out/Gene/filtered/matrix.mtx'),
        features=join(config['aligned_fastq10x_dir'], "{sample10x}",
                      'Solo.out/Gene/filtered/features.tsv'),
        barcodes=join(config['aligned_fastq10x_dir'], "{sample10x}",
                      'Solo.out/Gene/filtered/barcodes.tsv'),
        bam_alignments=join(config['aligned_fastq10x_dir'], "{sample10x}",
                            'Aligned.sortedByCoord.out.bam')
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
            '--outSAMattributes', 'NH', 'HI', 'nM', 'AS', 'MD', 'CB', 'UB',
            '--outSAMtype', 'BAM', 'SortedByCoordinate',
            '--runThreadN', str(threads),
            '--outFileNamePrefix', params.outdir,
            '--scoreGapNoncan', str(config['scoreGapNoncan']),
            '--scoreGapGCAG', str(config['scoreGapGCAG']),
            '--scoreGapATAC', str(config['scoreGapATAC']),
            ]
        print(f"Running STARsolo with following command:\n{' '.join(cmds)}")
        os.makedirs(params.outdir, exist_ok=True)
        subprocess.check_call(cmds)


rule get_cb_whitelist_10x:
    """Get whitelisted 10X cellbarcodes."""
    output: config['cb_whitelist_10x']
    params: url=config['cb_whitelist_10x_url']
    shell:
        """
        if [[ {params.url} == *.gz ]]
        then
            wget -O - {params.url} | gunzip -c > {output}
        else
            wget -O - {params.url} > {output}
        fi
        """
