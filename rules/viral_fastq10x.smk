"""Rules related viral reads in aligned 10x Illumina transcriptomics."""


rule viral_barcodes_in_transcripts:
    """Extract viral barcodes from 10x transcriptomic alignments."""
    input:
        bam=join(config['aligned_fastq10x_dir'], "{expt}",
                 'Aligned.sortedByCoord.out.bam'),
        bai=join(config['aligned_fastq10x_dir'], "{expt}",
                 'Aligned.sortedByCoord.out.bam.bai'),
        cell_barcodes=join(config['aligned_fastq10x_dir'], "{expt}",
                           'Solo.out/GeneFull/filtered/barcodes.tsv'),
        viral_bc_locs=join(config['viral_fastq10x_dir'], 'viral_bc_locs.csv'),
        notebook='notebooks/viral_barcodes_in_transcripts.py.ipynb'
    output:
        viral_bc_by_cell_umi_csv=join(config['viral_fastq10x_dir'],
                                      "{expt}_viral_bc_by_cell_umi.csv.gz"),
    log:
        notebook=join(config['viral_fastq10x_dir'],
                      "{expt}_viral_barcodes_in_transcripts.ipynb")
    notebook:
        '../notebooks/viral_barcodes_in_transcripts.py.ipynb'


rule viral_tags_in_transcripts:
    """Extract viral tags from 10x transcriptomic alignments."""
    input:
        bam=join(config['aligned_fastq10x_dir'], "{expt}",
                 'Aligned.sortedByCoord.out.bam'),
        bai=join(config['aligned_fastq10x_dir'], "{expt}",
                 'Aligned.sortedByCoord.out.bam.bai'),
        cell_barcodes=join(config['aligned_fastq10x_dir'], "{expt}",
                           'Solo.out/GeneFull/filtered/barcodes.tsv'),
        viral_tag_locs=join(config['viral_fastq10x_dir'], 'viral_tag_locs.csv'),
        viral_tag_identities=config['viral_tag_identities'],
        notebook='notebooks/viral_tags_in_transcripts.py.ipynb'
    output:
        viral_tag_by_cell_umi_csv=join(config['viral_fastq10x_dir'],
                                       "{expt}_viral_tag_by_cell_umi.csv.gz"),
    log:
        notebook=join(config['viral_fastq10x_dir'],
                      "{expt}_viral_tags_in_transcripts.ipynb")
    notebook:
        '../notebooks/viral_tags_in_transcripts.py.ipynb'


rule viral_bc_locs:
    """Locations of viral barcodes in 1-based indexing."""
    input:
        viral_genbank=config['viral_genbank']
    output:
        viral_bc_locs=join(config['viral_fastq10x_dir'], 'viral_bc_locs.csv'),
    run:
        viral_bc_tups = []
        for s in Bio.SeqIO.parse(input.viral_genbank, 'genbank'):
            for f in s.features:
                if f.type == 'viral_barcode':
                    viral_bc_tups.append((s.id,
                                          int(f.location.start) + 1,
                                          int(f.location.end)))
        pd.DataFrame.from_records(viral_bc_tups,
                                  columns=['gene', 'start', 'end']
                                  ).to_csv(output.viral_bc_locs, index=False)


rule viral_tag_locs:
    """Locations of viral tags in 1-based indexing."""
    input:
        viral_genbank=config['viral_genbank']
    output:
        viral_tag_locs=join(config['viral_fastq10x_dir'], 'viral_tag_locs.csv'),
    run:
        viral_tag_tups = []
        for s in Bio.SeqIO.parse(input.viral_genbank, 'genbank'):
            for f in s.features:
                if 'tag' in f.type:
                    viral_tag_tups.append((s.id,
                                           f.type,
                                           int(f.location.start) + 1,
                                           int(f.location.end)))
        pd.DataFrame.from_records(viral_tag_tups,
                                  columns=['gene', 'tag_name', 'start', 'end']
                                  ).to_csv(output.viral_tag_locs, index=False)
