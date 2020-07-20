"""``snakemake`` rule that loads required modules on the cluster-new partition."""

rule load_modules:
    output:
        'modules_loaded.txt'
    envmodules:
        'bcl2fastq/1.8.4'
    shell:
        """
        touch modules_loaded.txt
        """