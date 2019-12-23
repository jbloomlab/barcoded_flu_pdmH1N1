"""``snakemake`` rules calling viral tags / barcodes in 10X Illumina reads."""


rule viral_fastq10x_coverage:
    """Coverage of 10X Illumina reads on viral tags and barcodes."""
    input:
        nb='notebooks/viral_fastq10x_coverage.ipynb',
        fastq10x_bams=expand(join(config['aligned_fastq10x_dir'],
                                  "{sample10x}",
                                  'Aligned.sortedByCoord.out.bam'),
                             sample10x=samples_10x),
        fastq10x_bais=expand(join(config['aligned_fastq10x_dir'],
                                  "{sample10x}",
                                  'Aligned.sortedByCoord.out.bam.bai'),
                             sample10x=samples_10x),
        viral_genbank=config['viral_genbank']
    output:
        viraltag_locs=join(config['viral_fastq10x_dir'], 'viraltag_locs.csv'),
        viralbc_locs=join(config['viral_fastq10x_dir'], 'viralbc_locs.csv'),
        nb=join(config['viral_fastq10x_dir'], 'viral_fastq10x_coverage.ipynb'),
        nb_html=report(join(config['viral_fastq10x_dir'],
                            'viral_fastq10x_coverage.html'),
                       caption='../viral_fastq10x_coverage.rst',
                       category='Viral tags and barcodes in 10X data')
    run:
        run_nb_to_html(input_nb=input.nb,
                       output_nb=output.nb,
                       parameters={
                            'samples_10x': samples_10x,
                            'input_fastq10x_bams': input.fastq10x_bams,
                            'input_fastq10x_bais': input.fastq10x_bais,
                            'input_viral_genbank': input.viral_genbank,
                            'output_viraltag_locs': output.viraltag_locs,
                            'output_viralbc_locs': output.viralbc_locs,
                            },
                       )

rule index_bam:
    """Index a BAM file using `samtools`."""
    input: "{bamfile_base}.bam"
    output: "{bamfile_base}.bam.bai"
    shell:
        """samtools index {input} {output}"""
