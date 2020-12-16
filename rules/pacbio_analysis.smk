"""Rules related to analysis viral pacbio data."""

rule strand_exchange:
    """Quantify strand exchange."""
    input:
        viral_tags=config['viral_tag_identities'],
        mutation_df=join(config['pacbio_dir'],
                         "{expt}_mutations_by_ccs.csv.gz"),
        notebook='notebooks/strand_exchange.py.ipynb'
    output:
        plot_strand_exchange=report(join(config['pacbio_dir'],
                                    "{expt}_plot_strand_exchange.svg"),
                                    caption='../report/strand_exchange.rst',
                                    category="{expt}"),
        mutation_and_tag_identity_df=join(config['pacbio_dir'],
                                          "{expt}_mutations_str_ex.csv.gz")
    conda: '../environment.yml'
    log:
        notebook=join(config['log_dir'],
                      "strand_exchange_{expt}.ipynb")
    notebook:
        '../notebooks/strand_exchange.py.ipynb'

rule align_pacbio:
    """Summarize and aggregate all PacBio runs for an experiment."""
    input:
        amplicons=config['viral_amplicons'],
        features=config['amplicon_features'],
        ccs_fastq=lambda wc: [join(config['pacbio_dir'],
                                   f"{expt_pacbio_run}_ccs_fixterm.fastq.gz")
                              for expt_pacbio_run in
                              expts.expt_pacbio_runs(wc.expt)],
        ccs_report=lambda wc: [join(config['pacbio_dir'],
                                    f"{expt_pacbio_run}_report.txt")
                               for expt_pacbio_run in
                               expts.expt_pacbio_runs(wc.expt)],
        amplicon_to_reference_df=config['amplicon_to_reference'],
        notebook='notebooks/align_pacbio.py.ipynb'
    params:
        runs=lambda wc: expts.expt_pacbio_runs(wc.expt),
        alignment_stats=join(config['align_and_parse_dir'],
                             "{expt}")
    output:
        plot_amplicons=report(join(config['pacbio_dir'],
                              "{expt}_amplicons.svg"),
                              caption='../report/align_pacbio.rst',
                              category="{expt}"),
        plot_alignment=report(join(config['pacbio_dir'],
                              "{expt}_pacbio_passing_alignments.svg"),
                              caption='../report/align_pacbio.rst',
                              category="{expt}"),
        mutation_df=join(config['pacbio_dir'],
                         "{expt}_mutations_by_ccs.csv.gz"),
    threads: config['max_cpus'],
    conda: '../environment.yml'
    log:
        notebook=join(config['log_dir'],
                      "align_pacbio_{expt}.ipynb")
    notebook:
        '../notebooks/align_pacbio.py.ipynb'

rule pacbio_fix_termini:
    """Remove duplicate primer sequences at PacBio CCS termini."""
    input:
        fastq_in="{fastq_base}_ccs.fastq.gz",
        terminal_seqs=config['amplicon_linearization_primers']
    output:
        fastq_out="{fastq_base}_ccs_fixterm.fastq.gz"
    params:
        n_mismatches=1  # allowed mismatches in each copy of terminal seq
    conda: '../environment.yml'
    log:
        join(config['log_dir'], "pacbio_fix_termini_{fastq_base}.log")
    shell:
        """
        python scripts/remove_duplicate_terminal_seqs.py \
            {input.fastq_in} \
            {output.fastq_out} \
            {input.terminal_seqs} \
            --n_mismatches {params.n_mismatches} \
            &> {log}
        """

rule ccs_summaries:
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
        notebook='notebooks/ccs_summaries.py.ipynb'
    params:
        runs=lambda wc: expts.expt_pacbio_runs(wc.expt)

    output:
        summary=report(join(config['pacbio_dir'],
                            "{expt}_ccs_summaries.svg"),
                       caption='../report/ccs_summaries.rst',
                       category="{expt}")
    conda: '../environment.yml'
    threads: config['max_cpus']
    log:
        notebook=join(config['log_dir'],
                      "ccs_summaries_{expt}.ipynb")
    notebook:
        '../notebooks/ccs_summaries.py.ipynb'

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
