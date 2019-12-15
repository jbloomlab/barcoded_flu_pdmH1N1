"""``snakemake`` rules related to creating the ``STAR` reference genome."""


rule get_cell_genome:
    """Get the FASTA file for the cellular genome."""
    output: join(config['genome_dir'], 'cell_genome.fasta')
    params: ftp=config['cell_genome_ftp']
    shell: "wget -O - {params.ftp} | gunzip -c > {output}"


rule get_cell_gtf:
    """Get the GTF file for the cellular genome."""
    output: join(config['genome_dir'], 'cell_gtf.gtf')
    params: ftp=config['cell_gtf_ftp']
    shell: "wget -O - {params.ftp} | gunzip -c > {output}"


rule make_refgenome:
    """Build the STAR cellular + viral reference genome."""
    input:
        cell_genome=join(config['genome_dir'], 'cell_genome.fasta'),
        cell_gtf=join(config['genome_dir'], 'cell_gtf.gtf'),
        viral_genome=config['viral_genome'],
        viral_gtf=config['viral_gtf']
    output:
        concat_gtf=join(config['genome_dir'], 'cell_and_virus_gtf.gtf'),
        genomeDir=directory(config['refgenome'])
    threads: config['max_cpus']
    shell:
        """
        cat {input.cell_gtf} {input.viral_gtf} > {output.concat_gtf}
        mkdir -p {output.genomeDir}
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {output.genomeDir} \
             --genomeFastaFiles {input.cell_genome} {input.viral_genome} \
             --sjdbGTFfile {output.concat_gtf}
        """
