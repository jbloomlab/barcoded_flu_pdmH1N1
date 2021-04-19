"""Rules for making and QC-ing 10X transcriptomics FASTQ files."""


rule qc_fastq10x:
    """Quality control analysis of the 10x FASTQ files."""
    input:
        qc_stats=lambda wc: [join(config['fastq10x_dir'],
                                  f"{expt_run10x}_qc_stats.csv")
                             for expt_run10x
                             in expts.expt_transcriptomic_runs(wc.expt)],
        notebook='notebooks/qc_fastq10x.py.ipynb'
    output:
        qc_plot=report(join(config['fastq10x_dir'], "{expt}_qc_fastq10x.svg"),
                       caption='../report/qc_fastq10x.rst',
                       category="{expt}")
    log:
        notebook=join(config['log_dir'], "qc_fastq10x_{expt}.ipynb")
    conda: '../environment.yml'
    notebook:
        '../notebooks/qc_fastq10x.py.ipynb'


rule mkfastq10x:
    """Make 10x FASTQ files from BCL runs using `cellranger mkfastq`."""
    output:
        fastqR1=join(config['fastq10x_dir'], "{run10x}_all_R1.fastq.gz"),
        fastqR2=join(config['fastq10x_dir'], "{run10x}_all_R2.fastq.gz"),
        mkfastq10x_dir=directory(join(config['mkfastq10x_dir'], "{run10x}")),
        qc_stats=join(config['fastq10x_dir'], "{run10x}_qc_stats.csv"),
        csv=temp("_mkfastq_{run10x}.csv"),
        mro_file=temp("__{run10x}.mro"),
    params:
        index=lambda wc: expts.transcriptomic_index(wc.run10x),
        lane=lambda wc: expts.transcriptomic_lane(wc.run10x),
        index_sequencing=lambda wc: expts.transcriptomic_index_sequencing(wc.run10x),
        bcl_folder=lambda wc: expts.transcriptomic_bcl_folder(wc.run10x),
    threads:
        config['max_cpus']
    log:
        log=join(config['log_dir'], "mkfastq10x_{run10x}.log")
    conda: '../environment.yml'
    script:
        '../scripts/mkfastq10x.py'
