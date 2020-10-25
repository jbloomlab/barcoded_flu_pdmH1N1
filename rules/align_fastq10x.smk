"""Rules related aligned the 10X Illumina transcriptomics FASTQ reads."""


rule qc_transcript_alignments:
    """Quality control summary of 10x transcriptomics alignments."""
    input:
        summary=join(config['aligned_fastq10x_dir'], "{expt}",
                     'Solo.out/GeneFull/Summary.csv'),
        umi_per_cell=join(config['aligned_fastq10x_dir'], "{expt}",
                          'Solo.out/GeneFull/UMIperCellSorted.txt'),
        notebook='notebooks/qc_transcript_alignments.py.ipynb'
    output:
        qc_plot=report(join(config['aligned_fastq10x_dir'], "{expt}",
                            'qc_transcript_alignments.svg'),
                       caption='../report/qc_transcript_alignments.rst',
                       category="{expt}")
    log:
        notebook=join(config['aligned_fastq10x_dir'], "{expt}",
                      'qc_transcript_alignments.ipynb')
    conda: '../environment.yml'
    notebook:
        '../notebooks/qc_transcript_alignments.py.ipynb'


rule align_fastq10x:
    """Align 10x transcriptomics FASTQ reads with ``STARsolo``."""
    input:
        cb_whitelist_10x=config['cb_whitelist_10x'],
        refgenome=config['refgenome'],
        fastqR1=lambda wc: [join(config['fastq10x_dir'],
                                 f"{expt_run10x}_all_R1.fastq.gz")
                            for expt_run10x in 
                            expts.expt_transcriptomic_runs(wc.expt)],
        fastqR2=lambda wc: [join(config['fastq10x_dir'],
                                 f"{expt_run10x}_all_R2.fastq.gz")
                            for expt_run10x in 
                            expts.expt_transcriptomic_runs(wc.expt)],
    output:
        summary=join(config['aligned_fastq10x_dir'], "{expt}",
                     'Solo.out/GeneFull/Summary.csv'),
        umi_per_cell=join(config['aligned_fastq10x_dir'], "{expt}",
                          'Solo.out/GeneFull/UMIperCellSorted.txt'),
        matrix=join(config['aligned_fastq10x_dir'], "{expt}",
                    'Solo.out/GeneFull/filtered/matrix.mtx'),
        features=join(config['aligned_fastq10x_dir'], "{expt}",
                      'Solo.out/GeneFull/filtered/features.tsv'),
        barcodes=join(config['aligned_fastq10x_dir'], "{expt}",
                      'Solo.out/GeneFull/filtered/barcodes.tsv'),
        bam_alignments=join(config['aligned_fastq10x_dir'], "{expt}",
                            'Aligned.sortedByCoord.out.bam'),
    params:
        readFilesIn=lambda wc, input: (','.join(input.fastqR2) +
                                       ' ' +
                                       ','.join(input.fastqR1)
                                       ),
        expect_ncells=lambda wc: expts.expect_ncells(wc.expt),
        outdir=lambda wc, output: output.summary[: output.summary.index('Solo')],
        soloCBlen=config['cb_len_10x'],
        soloUMIstart=config['cb_len_10x'] + 1,
        soloUMIlen=config['umi_len_10x'],
        scoreGapNoncan=config['scoreGapNoncan'],
        scoreGapGCAG=config['scoreGapGCAG'],
        scoreGapATAC=config['scoreGapATAC'],
    threads: config['max_cpus']
    log: join(config['log_dir'], "align_fastq10x_{expt}.log")
    conda: '../environment.yml'
    shell:
        """STAR \
            --genomeDir {input.refgenome} \
            --readFilesIn {params.readFilesIn} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.outdir} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI nM AS MD CB UB \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist {input.cb_whitelist_10x} \
            --soloCBstart 1 \
            --soloCBlen {params.soloCBlen} \
            --soloUMIstart {params.soloUMIstart} \
            --soloUMIlen {params.soloUMIlen} \
            --soloBarcodeReadLength 0 \
            --soloCBmatchWLtype 1MM_multi_pseudocounts \
            --soloStrand Forward \
            --soloFeatures GeneFull \
            --soloUMIdedup 1MM_All \
            --soloUMIfiltering MultiGeneUMI \
            --soloOutFileNames Solo.out/ \
                               features.tsv \
                               barcodes.tsv \
                               matrix.mtx \
            --soloCellFilter CellRanger2.2 \
                             {params.expect_ncells} \
                             0.99 \
                             10 \
            --soloOutFormatFeaturesGeneField3 Gene Expression \
            --runThreadN {threads} \
            --scoreGapNoncan {params.scoreGapNoncan} \
            --scoreGapGCAG {params.scoreGapGCAG} \
            --scoreGapATAC {params.scoreGapATAC} \
        &> {log}
        """


rule get_cb_whitelist_10x:
    """Get whitelisted 10X cellbarcodes."""
    output: config['cb_whitelist_10x']
    params: url=config['cb_whitelist_10x_url']
    log: join(config['log_dir'], 'get_cb_whitelist_10x.log')
    conda: '../environment.yml'
    shell:
        """
        if [[ {params.url} == *.gz ]]
        then
            wget -O - {params.url} | gunzip -c > {output} 2> {log}
        else
            wget -O - {params.url} > {output} 2> {log}
        fi
        """
