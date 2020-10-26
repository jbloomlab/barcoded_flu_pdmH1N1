"""Rules related to analysis viral pacbio data."""

rule aggregate_expt_pacbio:
    """Aggregate all PacBio runs for an experiment."""
    input:
        ccs_reports=lambda wc: [join(config['pacbio_dir'],
                                     f"{expt_pacbio_run}_ccs.fastq.gz")
                                for expt_pacbio_run in
                                expts.expt_pacbio_runs(wc.expt)],
    output:
        summary=join(config['pacbio_dir'], "{expt}_summary.svg"),
    shell:
        # This rule should somehow aggregate all the summaries to make
        # a SVG plot that gives the stats, and create a concatenated
        # FASTQ with all the CCSs for the experiment. Right now it is just
        # a stand-in rule that doesn't do that, but just "touches" the
        # files to create them as empty. However, it should be possible
        # to easily create the summaries using an `alignparse.ccs.Summaries`
        # object.
        """
        touch {output.summary}
        """

rule build_ccs:
    """Run PacBio ``ccs`` program to build CCSs from subreads."""
    input:
     	subreads=lambda wc: expts.pacbio_subreads(wc.pacbio_run)
    output:
     	ccs_fastq=join(config['pacbio_dir'], "{pacbio_run}_ccs.fastq.gz"),
     	ccs_report=join(config['pacbio_dir'],"{pacbio_run}_report.txt"),
    params:
        ccs_min_length=config['ccs_min_length'],
        ccs_max_length=config['ccs_max_length'],
        ccs_min_rq=config['ccs_min_rq'],
    threads: config['max_cpus']
    shell:
	    """
        ccs \
            --report-file {output.ccs_report} \
            --num-threads {threads} \
            --min-length {params.ccs_min_length} \
            --max-length {params.ccs_max_length} \
            --min-rq {params.ccs_min_rq} \
            {input.subreads} \
            {output.ccs_fastq}
		"""

