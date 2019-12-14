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

# extract Illumina 10X runs from config
illumina_runs_10x = (
    pd.read_csv(config['illumina_runs_10x'], comment='#')
    .assign(run10x=lambda x: x['sample'] + '-' + x['seqrun'].astype(str))
    .set_index('run10x')
    )
assert len(illumina_runs_10x) == illumina_runs_10x.index.nunique()


# Target rules ---------------------------------------------------------------

rule all:
    input:
        expand(join(config['fastq10x_dir'], "{run10x}_all_{read}.fastq.gz"),
               run10x=illumina_runs_10x.index, read=['R1', 'R2'])


# Load rules -----------------------------------------------------------------

include: 'rules/fastq10x.smk'
