"""Rules related to analysis viral pacbio data."""

rule build_ccs:
    """Run PacBio ``ccs`` program to build CCSs from subreads."""
    
    input:
     	subreads=lambda wc: expts.pacbio_subreads(wc.pacbio_run)

    output:
     	ccs_fastq = join(config['pacbio_dir'], "{pacbio_run}_ccs.fastq.gz"),
     	ccs_report = join(config['pacbio_dir'],"{pacbio_run}_report.txt"),


    threads: config['max_cpus']

    shell:
	    """
	     	ccs \
	     		--report-file {output.ccs_report} \
	     		--num-threads {threads} \
	     		{input.subreads} \
	     		{output.ccs_fastq}
		"""

