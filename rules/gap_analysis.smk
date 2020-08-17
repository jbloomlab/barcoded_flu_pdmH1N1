"""``snakemake`` rules related to analysis of gapped reads."""


rule analyze_gaps:
    """Analyze gaps in 10X sequencing reads."""
    input:
        nb='notebooks/gap_analysis.ipynb',
        fastq10x_bams=expand(join(config['aligned_fastq10x_dir'],
                                  "{sample10x}/Aligned.sortedByCoord.out.bam"),
                             sample10x=samples_10x),
        fastq10x_bais=expand(join(config['aligned_fastq10x_dir'],
                                  "{sample10x}",
                                  'Aligned.sortedByCoord.out.bam.bai'),
                             sample10x=samples_10x),
        viral_genbank=config['viral_genbank'],
        viraltag_locs=join(config['viral_fastq10x_dir'], 'viraltag_locs.csv'),
        viralbc_locs=join(config['viral_fastq10x_dir'], 'viralbc_locs.csv')
    output:
        nb=join(config['viral_fastq10x_dir'], 'gap_analysis.ipynb'),
        nb_html=report(join(config['viral_fastq10x_dir'],
                            'gap_analysis.html'),
                       caption='../report/gap_analysis.rst',
                       category='Viral tags and barcodes in 10X data'),
        gapped_reads=join(config['viral_fastq10x_dir'], 'gapped_reads.csv'),
        gapped_reads_summary=join(config['viral_fastq10x_dir'], 'gapped_reads_summary.csv')
    run:
        run_nb_to_html(input_nb=input.nb,
                   output_nb=output.nb,
                   parameters={
                        'samples_10x': samples_10x,
                        'input_fastq10x_bams': input.fastq10x_bams,
                        'input_fastq10x_bais': input.fastq10x_bais,
                        'input_viral_genbank': input.viral_genbank,
                        'input_viraltag_locs': input.viraltag_locs,
                        'input_viralbc_locs': input.viralbc_locs,
                        'output_gapped_reads': output.gapped_reads,
                        'output_gapped_reads_summmary': output.gapped_reads_summary,
                        },
                   )
