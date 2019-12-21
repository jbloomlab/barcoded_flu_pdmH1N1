"""``snakemake`` rules for calling viral barcodes from 10X Illumina reads."""


#rule get_viralbc_fastq10x:
#    """Get 10X reads aligned to barcoded viral gene with valid cell barcode."""
#    input:
#        cell_barcodes=join(config['aligned_fastq10x_dir'], "{sample10x}",
#                           'Solo.out/Gene/filtered/barcodes.tsv'),
#        bam_alignments=join(config['aligned_fastq10x_dir'], "{sample10x}",
#                            'Aligned.sortedByCoord.out.bam')
#    params:
#        bc_viral_gene="{bc_viral_gene}"
#    output:
#        fasta=join(config['viralbc_fastq10x_dir'],
#                   "{sample10x}_{bc_viral_gene}_reads.fasta.gz")
#    shell:
        # see here to understand shell command:
        # https://kb.10xgenomics.com/hc/en-us/articles/360022448251-Is-there-way-to-filter-the-BAM-file-produced-by-10x-pipelines-with-a-list-of-barcodes-
#        """
#        samtools view {input.bam_alignments} "
#        """

rule index_bam:
    """Index a BAM file using `samtools`."""
    input: "{bamfile_base}.bam"
    output: "{bamfile_base}.bam.bai"
    shell:
        """samtools index {input} {output}"
