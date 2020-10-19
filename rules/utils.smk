"""General utility rules."""


rule index_bam:
    """Index a BAM file using `samtools`."""
    input: "{bamfile_base}.bam"
    output: "{bamfile_base}.bam.bai"
    conda: '../environment.yml'
    shell:
        """samtools index {input} {output}"""
