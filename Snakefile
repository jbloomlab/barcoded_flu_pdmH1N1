"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------

import glob
import json
import os
from os.path import join, basename
import re
import shutil
import subprocess

import pandas as pd

from pymodules.jupnb import run_nb_to_html


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

# list of all 10X samples
samples_10x = illumina_runs_10x['sample'].unique().tolist()


# Target rules ---------------------------------------------------------------

rule all:
    input:
        join(config['fastq10x_dir'], 'fastq10x_qc_analysis.html'),
        join(config['aligned_fastq10x_dir'], 'align_fastq10x_summary.html'),
        expand(join(config['analysis_dir'], 
                    "{sample10x}_analyze_cell_gene_matrix.html"),
               sample10x=samples_10x)


# Set up report  -------------------------------------------------------------

report: 'report/workflow.rst'


# Load rules -----------------------------------------------------------------

include: 'rules/analysis.smk'
include: 'rules/align_fastq10x.smk'
include: 'rules/star_refgenome.smk'
include: 'rules/fastq10x.smk'
