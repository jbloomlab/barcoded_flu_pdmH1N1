"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import json
from os.path import join, basename
import re
import shutil
import subprocess

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: all

# Global variables extracted from config --------------------------------------
illumina_runs_10x = (
    pd.read_csv(config['illumina_runs_10x'], comment='#')
    .assign(run10x=lambda x: x['sample'] + '-' + x['seqrun'].astype(str))
    .set_index('run10x')
    )
assert len(illumina_runs_10x) == illumina_runs_10x.index.nunique()

# Global variables defined statically ----------------------------------------
genome=['cell_genome',
        'spikein_genome']

# Rules -----------------------------------------------------------------------

rule all:
    input:
        expand(join(config['fastq10x_dir'], "{run10x}_all_{read}.fastq.gz"),
               run10x=illumina_runs_10x.index, read=['R1', 'R2'])

ruleorder: filter_gtf > get_gtf
        
rule get_genome:
    """Download genome assembly in a FASTA file."""
    input:
    output:
        join(config['genome_dir'], "{genome}.fasta")
    params:
        ftp=lambda wildcards: config[wildcards.genome]['fasta']
    shell:
        "wget -O - {params.ftp} | gunzip -c > {output}"
        
rule get_gtf:
    """Download gene sets in a GTF file."""
    input:
    output:
        join(config['genome_dir'], "{genome}.gtf")
    params:
        ftp=lambda wildcards: config[wildcards.genome]['gtf']
    shell:
        "wget -O - {params.ftp} | gunzip -c > {output}"

rule filter_gtf:
    """Filter GTF files to only features of interest (e.g. protein-coding genes)."""
    input:
        join(config['genome_dir'], "{genome}.gtf")
    output:
        join(config['genome_dir'], "{genome}_filtered.gtf")
    shell:
        """
        cellranger mkgtf {input} {output} \
            --attribute=gene_biotype:protein_coding
        """


rule make_fastq10x:
    """Make 10X FASTQ files from BCL runs using `cellranger mkfastq`."""
    output:
        fastqR1=join(config['fastq10x_dir'], "{run10x}_all_R1.fastq.gz"),
        fastqR2=join(config['fastq10x_dir'], "{run10x}_all_R2.fastq.gz"),
        mkfastq10x_dir=directory(join(config['mkfastq10x_dir'], "{run10x}")),
        qc_stats=join(config['fastq10x_dir'], "{run10x}_qc_stats.csv"),
        csv=temp("_mkfastq_{run10x}.csv"),  # input for cellranger mkfastq
        mro=temp("__{run10x}.mro"),  # created by cellranger mkfastq
    params:
        run10x="{run10x}"
    threads:
        config['max_cpus']
    run:
        # write CSV file for `cellranger mkfastq`
        with open(output.csv, 'w') as f:
            f.write('Lane,Sample,Index\n' +
                    ','.join(map(str,
                             [illumina_runs_10x.at[params.run10x, 'lane'],
                              params.run10x,
                              illumina_runs_10x.at[params.run10x, 'index']])))

        # run `cellranger mkfastq`
        cmds = ['cellranger', 'mkfastq',
                '--run', illumina_runs_10x.at[params.run10x, 'bcl_folder'],
                '--id', params.run10x,  # output directory name
                '--csv', output.csv,
                '--delete-undetermined',
                '--qc',
                f"--localcores={threads}",
                ]
        print(f"\nRunning the following commands:\n{' '.join(cmds)}\n")
        subprocess.check_call(cmds)

        # move `cellranger mkfastq` output to desired location
        print(f"\nMoving `cellranger mkfastq` output from {params.run10x} "
              f"to {output.mkfastq10x_dir}\n")
        shutil.move(params.run10x, output.mkfastq10x_dir)

        # get names of R1 and R2 FASTQ files from `cellranger mkfastq` output
        fastq_glob = f"{output.mkfastq10x_dir}/outs/fastq_path/*/*/*.fastq.gz"
        fastqs = sorted(glob.glob(fastq_glob))
        fastqregex = re.compile(r'_(?P<read>R1|R2|I1)_\d{3}\.fastq\.gz')
        r1s = [f for f in fastqs if fastqregex.search(f).group('read') == 'R1']
        r2s = [f for f in fastqs if fastqregex.search(f).group('read') == 'R2']
        print(f"\nFASTQ files matching {fastq_glob}:\n"
              f"  R1 files: {','.join(map(basename, r1s))}\n"
              f"  R2 files: {','.join(map(basename, r2s))}\n")
        assert len(r1s) == len(r2s)

        # concatenate all R1 and R2 files into merged FASTQs for run
        print(f"\nCreating merged FASTQs {output.fastqR1}, {output.fastqR2}\n")
        for outfq, fqlist in [(output.fastqR1, r1s), (output.fastqR2, r2s)]:
            with open(outfq, 'wb') as f:
                subprocess.call(['cat'] + fqlist, stdout=f)

        # extract sample QC stats from `cellranger mkfastq` JSON into CSV
        qc_json = join(config['mkfastq10x_dir'], params.run10x,
                       'outs/qc_summary.json')
        print(f"\nExtracting QC stats from {qc_json} to {output.qc_stats}\n")
        with open(qc_json) as f:
            (pd.Series(json.load(f)['sample_qc'][params.run10x]['all'])
             .to_csv(output.qc_stats, header=False)
             )
