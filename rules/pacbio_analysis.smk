"""Rules related to analysis viral pacbio data."""

rule build_ccs:
    """Run PacBio ``ccs`` program to build CCSs from subreads."""
    
    input:
     	subreads=lambda wc: expts.pacbio_subreads(wc.expt)

    output:
     	ccs_fastq = join(config['pacbio_dir'], "{expt}_ccs.fastq.gz"),
     	ccs_report = join(config['pacbio_dir'],"{expt}_report.txt"),

    threads: config['max_cpus']

    shell:
	    """
	     	ccs \
	     		--report-file {output.ccs_report}
	     		--num-threads {threads}
	     		{params.subreads} \
	     		{output.ccs_fastq}
		"""

