"""``snakemake`` rules calling viral tags / barcodes in 10X Illumina reads."""


rule count_viralbc_fastq10x:
    """Count viral barcodes for each gene / cell from the 10X Illumina data."""
    input:
        nb='notebooks/count_viralbc_fastq10x.ipynb',
        fastq10x_bam=join(config['aligned_fastq10x_dir'],
                          "{sample10x}/Aligned.sortedByCoord.out.bam"),
        fastq10x_bai=join(config['aligned_fastq10x_dir'],
                          "{sample10x}/Aligned.sortedByCoord.out.bam.bai"),
        viralbc_locs=join(config['viral_fastq10x_dir'], 'viralbc_locs.csv'),
        cellbarcodes=join(config['aligned_fastq10x_dir'], "{sample10x}"
                          'Solo.out/Gene/filtered/barcodes.tsv')
    output:
        nb=join(config['viral_fastq10x_dir'],
                "count_viralbc_fastq10x-{sample10x}.ipynb"),
        nb_html=report(join(config['viral_fastq10x_dir'],
                            "count_viralbc_fastq10x-{sample10x}.html"),
                       caption='../report/count_viralbc_fastq10x.rst',
                       category='Viral tags and barcodes in 10X data'),
        viralbc_counts=join(config['viral_fastq10x_dir'],
                            "viralbc_counts_{sample10x}.csv")
    run:
        run_nb_to_html(
                input_nb=input.nb,
                output_nb=output.nb,
                parameters={
                    'input_fastq10x_bam': input.fastq10x_bam,
                    'input_fastq10x_bai': input.fastq10x_bai,
                    'input_viralbc_locs': input.viraltag_locs,
                    'input_cellbarcodes': input.cellbarcodes,
                    'output_viralbc_counts': output.viraltag_counts,
                    },
                )


rule count_viraltags_fastq10x:
    """Count viral tags for each gene and cell from the 10X Illumina data."""
    input:
        nb='notebooks/count_viraltags_fastq10x.ipynb',
        fastq10x_bam=join(config['aligned_fastq10x_dir'],
                          "{sample10x}/Aligned.sortedByCoord.out.bam"),
        fastq10x_bai=join(config['aligned_fastq10x_dir'],
                          "{sample10x}/Aligned.sortedByCoord.out.bam.bai"),
        viraltag_locs=join(config['viral_fastq10x_dir'], 'viraltag_locs.csv'),
        viraltag_identities=config['viraltag_identities'],
        cellbarcodes=join(config['aligned_fastq10x_dir'], "{sample10x}",
                          'Solo.out/Gene/filtered/barcodes.tsv')
    output:
        nb=join(config['viral_fastq10x_dir'],
                "count_viraltags_fastq10x-{sample10x}.ipynb"),
        nb_html=report(join(config['viral_fastq10x_dir'],
                            "count_viraltags_fastq10x-{sample10x}.html"),
                       caption='../report/count_viraltags_fastq10x.rst',
                       category='Viral tags and barcodes in 10X data'),
        viraltag_counts=join(config['viral_fastq10x_dir'],
                             "viraltag_counts_{sample10x}.csv")
    run:
        run_nb_to_html(
                input_nb=input.nb,
                output_nb=output.nb,
                parameters={
                    'input_fastq10x_bam': input.fastq10x_bam,
                    'input_fastq10x_bai': input.fastq10x_bai,
                    'input_viraltag_locs': input.viraltag_locs,
                    'input_viraltag_identities': input.viraltag_identities,
                    'input_cellbarcodes': input.cellbarcodes,
                    'output_viraltag_counts': output.viraltag_counts,
                    },
                )


rule viral_fastq10x_coverage:
    """Coverage of 10X Illumina reads on viral tags and barcodes."""
    input:
        nb='notebooks/viral_fastq10x_coverage.ipynb',
        fastq10x_bams=expand(join(config['aligned_fastq10x_dir'],
                                  "{sample10x}/Aligned.sortedByCoord.out.bam"),
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
                       caption='../report/viral_fastq10x_coverage.rst',
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
