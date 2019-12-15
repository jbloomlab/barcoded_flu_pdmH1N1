"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------

import glob
import json
from os.path import join, basename
import re
import shutil
import subprocess

# see here: https://stackoverflow.com/a/29172195
import matplotlib
matplotlib.use('Agg')

import mizani

import pandas as pd

from plotnine import *


# Configuration  --------------------------------------------------------------

configfile: 'config.yaml'

# run "quick" rules locally:
localrules: all,
            fastq10x_qc_stats,

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
               run10x=illumina_runs_10x.index, read=['R1', 'R2']),
        join(config['fastq10x_dir'], 'fastq10x_qc_stats.svg'),
        join(config['fastq10x_dir'], 'fastq10x_qc_stats.csv'),
        config['refgenome'],
        config['cb_whitelist_10x']


# Set up report  -------------------------------------------------------------

report: 'report/workflow.rst'


# Load rules -----------------------------------------------------------------

include: 'rules/align_fastq10x.smk'
include: 'rules/star_refgenome.smk'
include: 'rules/fastq10x.smk'
