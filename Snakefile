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

import pymodules.experiments


# Configuration  --------------------------------------------------------------

configfile: 'config.yaml'

expts = pymodules.experiments.Experiments(config['experiments'])

localrules: all


# Target rules ---------------------------------------------------------------

rule all:
    input:
#        expand(join(config['mkfastq10x_dir'],
#                    "{transcriptomic_run}_qc_stats.csv",
#               transcriptomic_run=experiments_config.transcriptomic_runs)
        expand(join(config['fastq10x_dir'], "{expt}_fastq10x_qc.svg"),
               expt=expts.experiments),
#        join(config['aligned_fastq10x_dir'], 'align_fastq10x_summary.html'),
#        join(config['aligned_fastq10x_dir'],
#             'fastq10x_transcript_coverage.html'),
#        join(config['viral_fastq10x_dir'], 'viral_fastq10x_coverage.html'),
#        join(config['viral_fastq10x_dir'], 'gap_analysis.html'),
#        expand(join(config['viral_fastq10x_dir'],
#                    "count_viraltags_fastq10x-{sample10x}.html"),
#               sample10x=experiments_config.experiments),
#        expand(join(config['viral_fastq10x_dir'],
#                    "count_viralbc_fastq10x-{sample10x}.html"),
#               sample10x=experiments_config.experiments),
#        expand(join(config['analysis_dir'], 
#                    "{sample10x}_analyze_cell_gene_matrix.html"),
#               sample10x=experiments_config.experiments),


# Set up report  -------------------------------------------------------------

report: 'report/workflow.rst'


# Load rules -----------------------------------------------------------------

#include: 'rules/analysis.smk'
#include: 'rules/viral_fastq10x.smk'
#include: 'rules/align_fastq10x.smk'
#include: 'rules/star_refgenome.smk'
include: 'rules/fastq10x.smk'
#include: 'rules/gap_analysis.smk'
