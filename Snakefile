"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------

import glob
import json
import os
from os.path import join, basename
import re
import shutil
import subprocess

import Bio.SeqIO

import pandas as pd

from pymodules.jupnb import run_nb_to_html


# Configuration  --------------------------------------------------------------

configfile: 'config.yaml'

# run "quick" rules locally:
localrules: all


# Global variables used in analysis---------------------------------------------

# Get Illumina 10X runs and samples
illumina_runs_10x = (
    pd.read_csv(config['illumina_runs_10x'], comment='#')
    .assign(run10x=lambda x: x['sample'] + '-' + x['seqrun'].astype(str))
    .set_index('run10x')
    )
assert len(illumina_runs_10x) == illumina_runs_10x.index.nunique()
samples_10x = illumina_runs_10x['sample'].unique().tolist()  # list 10X samples

# list of barcoded viral genes, extracted from viral genome Genbank file
bc_viral_genes = []
for seqrecord in Bio.SeqIO.parse(config['viral_genome'], 'genbank'):
    nbcs = len([f for f in seqrecord.features if f.type == 'viral_barcode'])
    if nbcs == 1:
        bc_viral_genes.append(seqrecord.id)
    elif nbcs > 0:
        raise ValueError(f"multiple viral barcodes for {seqrecord.id}")


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
include: 'rules/viral_fastq10x.smk'
include: 'rules/align_fastq10x.smk'
include: 'rules/star_refgenome.smk'
include: 'rules/fastq10x.smk'
