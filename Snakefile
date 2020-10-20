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

import yaml


# Configuration  --------------------------------------------------------------

configfile: 'config.yaml'

expts = pymodules.experiments.Experiments(config['experiments'])

# get possible viral tag identities
with open(config['viral_tag_identities']) as f:
    viral_tags = sorted({tag_variant for gene_tags in yaml.safe_load(f).values()
                         for tags in gene_tags.values()
                         for tag_variant in tags})

# get viral genes
viral_genes = [s.id for s in Bio.SeqIO.parse(config['viral_genbank'],
                                             'genbank')]
assert viral_genes == [s.id for s in Bio.SeqIO.parse(config['viral_genome'],
                                                     'fasta')]

# get barcoded viral genes
barcoded_viral_genes = [s.id for s in Bio.SeqIO.parse(config['viral_genbank'],
                                                      'genbank')
                        if any(f.type == 'viral_barcode' for f in s.features)]


# Target rules ---------------------------------------------------------------

localrules: all

rule all:
    input:
        expand(join(config['fastq10x_dir'], "{expt}_qc_fastq10x.svg"),
               expt=expts.experiments),
        expand(join(config['aligned_fastq10x_dir'], "{expt}",
                    'qc_transcript_alignments.svg'),
               expt=expts.experiments),
        expand(join(config['viral_fastq10x_dir'],
                    "{expt}_viral_tag_by_cell.svg"),
               expt=expts.experiments),
        expand(join(config['viral_fastq10x_dir'],
                    "{expt}_viral_bc_by_cell.svg"),
               expt=expts.experiments),
        expand(join(config['viral_fastq10x_dir'],
                    "{expt}_viral_transcript_coverage.svg"),
               expt=expts.experiments),
#        join(config['viral_fastq10x_dir'], 'viral_fastq10x_coverage.html'),
#        join(config['viral_fastq10x_dir'], 'viral_fastq10x_coverage.html'),
#        join(config['viral_fastq10x_dir'], 'gap_analysis.html'),
#        expand(join(config['viral_fastq10x_dir'],
#                    "count_viraltags_fastq10x-{sample10x}.html"),
#               sample10x=experiments_config.experiments),
#        expand(join(config['viral_fastq10x_dir'],
#                    "count_viralbc_fastq10x-{sample10x}.html"),
#               sample10x=experiments_config.experiments),


# Set up report  -------------------------------------------------------------

report: 'report/workflow.rst'


# Load rules -----------------------------------------------------------------

include: 'rules/viral_fastq10x.smk'
include: 'rules/align_fastq10x.smk'
include: 'rules/star_refgenome.smk'
include: 'rules/fastq10x.smk'
include: 'rules/utils.smk'
