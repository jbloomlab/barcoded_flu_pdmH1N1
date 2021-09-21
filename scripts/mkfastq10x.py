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

# build samplesheet from available options
data = ['[Data],,,,,,,,\n',]
values = []

data.extend(['Lane,'])
values.extend([f'{snakemake.params.lane},'.strip("'")])

data.extend(['Sample_ID,','Sample_Name,','Sample_Plate,','Sample_Well,'])
values.extend([f'{snakemake.wildcards.run10x},',
               f'{snakemake.wildcards.run10x}',
               ',,,'])

data.extend(['I7_Index_ID,','index,'])
if snakemake.params.index_sequencing == 'none':
    values.extend([',,'])
else:
    values.extend([f'{snakemake.params.index},',f'{snakemake.params.index},'])
if snakemake.params.index_sequencing == 'dual':
    data.extend(['I5_Index_ID,','index2,'])
    values.extend([f'{snakemake.params.index},',f'{snakemake.params.index},'])

data.extend(['Sample_Project,','Description\n'])
values.extend(['project,,'])

# convert samplesheet to string
samplesheet_contents = []
for part in [data, values]:
    samplesheet_contents.append(''.join(part))
samplesheet_contents = ''.join(samplesheet_contents)

print(f"Writing samplesheet to {snakemake.output.samplesheet}")
with open(snakemake.output.samplesheet, 'w') as f:
    f.write(samplesheet_contents)

# run `cellranger mkfastq`
cmds = ['cellranger', 'mkfastq',
        '--run', snakemake.params.bcl_folder,
        '--id', snakemake.wildcards.run10x,  # output directory name
        '--samplesheet', snakemake.output.samplesheet,
        '--delete-undetermined',
        '--qc',
        f"--localcores={snakemake.threads}",
        ]
if snakemake.params.index_sequencing == 'none':  # must use --lanes flag if no sample index
    cmds.extend(['--lanes'])
    if snakemake.params.lane == '*':
        cmds.extend(['1,2'])
    else:
        cmds.extend([str(snakemake.params.lane)])
if snakemake.params.index_sequencing == 'single':
    cmds.extend(['--force-single-index'])
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
                          'outs/fastq_path/**/*.fastq.gz')
fastqs = sorted(glob.glob(fastq_glob, recursive=True))
fastqregex = re.compile(r'_(?P<read>R1|R2|I1|I2)_\d{3}\.fastq\.gz')
r1s = [f for f in fastqs if fastqregex.search(f).group('read') == 'R1']
r2s = [f for f in fastqs if fastqregex.search(f).group('read') == 'R2']
print(f"\nFASTQ files matching {fastq_glob}:\n"
      f"  R1 files: {','.join(map(os.path.basename, r1s))}\n"
      f"  R2 files: {','.join(map(os.path.basename, r2s))}\n")
assert len(r1s) == len(r2s) > 0

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
