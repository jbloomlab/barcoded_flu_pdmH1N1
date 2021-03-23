"""Make 10x FASTQ files from BCL runs using `cellranger mkfastq`."""


import glob
import json
import os
import re
import shutil
import subprocess
import sys

import pandas as pd


print(f"Writing output and errors to {snakemake.log.log}")
f = open(snakemake.log.log, 'w')
sys.stdout = f
sys.stderr = f

print(f"Writing CSV to {snakemake.output.csv}")
with open(snakemake.output.csv, 'w') as f:
    f.write('Lane,Sample,Index\n' + ','.join([snakemake.params.lane,
                                              snakemake.wildcards.run10x,
                                              snakemake.params.index,
                                              ])
            )

# run `cellranger mkfastq`
cmds = ['cellranger', 'mkfastq',
        '--run', snakemake.params.bcl_folder,
        '--id', snakemake.wildcards.run10x,  # output directory name
        '--csv', snakemake.output.csv,
        '--delete-undetermined',
        '--qc',
        f"--localcores={snakemake.threads}",
        ]
if 'GA' in snakemake.params.index:
    cmds.extend(['--use-bases-mask', 'Y*,I8n*,Y*'])
print(f"\nRunning the following commands:\n{' '.join(cmds)}\n")
subprocess.check_call(cmds)

# move `cellranger mkfastq` output to desired location
print(f"\nMoving `cellranger mkfastq` output from {snakemake.wildcards.run10x}"
      f" to {snakemake.output.mkfastq10x_dir}\n")
shutil.move(snakemake.wildcards.run10x, snakemake.output.mkfastq10x_dir)

# get names of R1 and R2 FASTQ files from `cellranger mkfastq` output
fastq_glob = os.path.join(snakemake.output.mkfastq10x_dir,
                          'outs/fastq_path/*/*/*.fastq.gz')
fastqs = sorted(glob.glob(fastq_glob))
fastqregex = re.compile(r'_(?P<read>R1|R2|I1)_\d{3}\.fastq\.gz')
r1s = [f for f in fastqs if fastqregex.search(f).group('read') == 'R1']
r2s = [f for f in fastqs if fastqregex.search(f).group('read') == 'R2']
print(f"\nFASTQ files matching {fastq_glob}:\n"
      f"  R1 files: {','.join(map(os.path.basename, r1s))}\n"
      f"  R2 files: {','.join(map(os.path.basename, r2s))}\n")
assert len(r1s) == len(r2s) == len(fastqs) / 3 > 0

# concatenate all R1 and R2 files into merged FASTQs for run
print(f"\nCreating merged FASTQs {snakemake.output.fastqR1}, "
      f"{snakemake.output.fastqR2}\n")
for outfq, fqlist in [(snakemake.output.fastqR1, r1s),
                      (snakemake.output.fastqR2, r2s)]:
    with open(outfq, 'wb') as f:
        subprocess.call(['cat'] + fqlist, stdout=f)

# extract sample QC stats from `cellranger mkfastq` JSON into CSV
qc_json = os.path.join(snakemake.output.mkfastq10x_dir, 'outs/qc_summary.json')
print(f"\nExtracting QC stats from {qc_json} to {snakemake.output.qc_stats}\n")
with open(qc_json) as f:
    (pd.Series(json.load(f)['sample_qc'][snakemake.wildcards.run10x]['all'])
     .to_csv(snakemake.output.qc_stats, header=False)
     )
