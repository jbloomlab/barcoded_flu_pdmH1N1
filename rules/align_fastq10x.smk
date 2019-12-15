"""``snakemake`` rules related aligned the Illumina FASTQ 10X reads."""


rule get_cb_whitelist_10x:
    """Get whitelisted 10X cellbarcodes."""
    output: config['cb_whitelist_10x']
    params: ftp=config['cb_whitelist_10x_ftp']
    shell: "wget -O - {params.ftp} | gunzip -c > {output}"
