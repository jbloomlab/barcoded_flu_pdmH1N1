"""Rules related to analysis viral pacbio data."""

rule ccs_summarize:
    """Summarize and aggregate all PacBio runs for an experiment."""
    input:
        ccs_fastq=lambda wc: [join(config['pacbio_dir'],
                                   f"{expt_pacbio_run}_ccs.fastq.gz")
                              for expt_pacbio_run in
                              expts.expt_pacbio_runs(wc.expt)],

        ccs_report=lambda wc: [join(config['pacbio_dir'],
                                    f"{expt_pacbio_run}_report.txt")
                               for expt_pacbio_run in
                               expts.expt_pacbio_runs(wc.expt)],
        notebook='notebooks/summarize_pacbio.py.ipynb'
    params:
        runs=lambda wc: expts.expt_pacbio_runs(wc.expt)

    output:
        summary=join(config['pacbio_dir'], "{expt}_ccs_summary.svg"),
    conda: '../environment.yml'
    log:
        notebook=join(config['log_dir'],
                      "{expt}_summarize_pacbio.ipynb")
    notebook:
        '../notebooks/summarize_pacbio.py.ipynb'

rule build_ccs:
    """Run PacBio ``ccs`` program to build CCSs from subreads."""
    input:
        subreads=lambda wc: expts.pacbio_subreads(wc.pacbio_run)
    output:
        ccs_fastq=join(config['pacbio_dir'], "{pacbio_run}_ccs.fastq.gz"),
        ccs_report=join(config['pacbio_dir'], "{pacbio_run}_report.txt"),
    params:
        ccs_min_length=config['ccs_min_length'],
        ccs_max_length=config['ccs_max_length'],
        ccs_min_rq=config['ccs_min_rq'],
    threads: config['max_cpus'],
    conda: '../environment.yml'
    log: join(config['log_dir'], "build_ccs_{pacbio_run}.log")
    shell:
        """
        ccs \
            --report-file {output.ccs_report} \
            --log-level INFO \
            --log-file {log} \
            --num-threads {threads} \
            --min-length {params.ccs_min_length} \
            --max-length {params.ccs_max_length} \
            --min-rq {params.ccs_min_rq} \
            {input.subreads} \
            {output.ccs_fastq}
        """
