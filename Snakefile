"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import os
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
    .assign(
        run10x=lambda x: x['sample'] + '-' + x['seqrun'].astype(str),
        )
    .set_index('run10x')
    )
assert len(illumina_runs_10x) == illumina_runs_10x.index.nunique()

# Rules -----------------------------------------------------------------------
rule:
    input:
        expand(os.path.join(config['fastq10x_dir'], "{run10x}_{r}.fastq.gz"),
               run10x=illumina_runs_10x.index, r=['R1', 'R2'])

rule mkfastq:
    """10X FASTQ files from BCL runs using `cellranger mkfastq`."""
    output:
        fastqR1=os.path.join(config['fastq10x_dir'], "{run10x}_R1.fastq.gz"),
        fastqR2=os.path.join(config['fastq10x_dir'], "{run10x}_R2.fastq.gz"),
        mkfastq10x_dir=directory(os.path.join(config['mkfastq10x_dir'],
                                              "{run10x}")),
        csv=temp("{run10x}.csv"),  # used as input by cellranger mkfastq
        mro=temp("__{run10x}.mro"),  # created by cellranger mkfastq
    params:
        run10x="{run10x}"
    run:
        # write CSV file for `cellranger mkfastq`
        with open(output.csv, 'w') as f:
            f.write('Lane,Sample,Index\n' +
                    ','.join(map(str,
                             [illumina_runs_10x.at[params.run10x, 'lane'],
                              illumina_runs_10x.at[params.run10x, 'sample'],
                              illumina_runs_10x.at[params.run10x, 'index']])))
        # run `cellranger mkfastq`
        cmds = ['cellranger', 'mkfastq',
                '--run', illumina_runs_10x.at[params.run10x, 'bcl_folder'],
                '--id', params.run10x,  # output directory name
                '--csv', output.csv,
                '--delete-undetermined',
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
              f"  R1 files: {','.join(map(os.path.basename, r1s))}\n"
              f"  R2 files: {','.join(map(os.path.basename, r2s))}\n")
        assert len(r1s) == len(r2s)
        # concatenate all R1 and R2 files into merged FASTQ for run
        print(f"\nCreating merged FASTQs {output.fastqR1} and {output.fastqR2}\n")
        for outfq, fqlist in [(output.fastqR1, r1s), (output.fastqR2, r2s)]:
            with open(outfq, 'wb') as f:
                subprocess.call(['cat'] + fqlist, stdout=f)
