"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------

import os
from os.path import join

import Bio.SeqIO

import pymodules.experiments

import yaml


# Configuration  --------------------------------------------------------------

configfile: 'config.yaml'

# get possible viral tag identities
with open(config['viral_tag_identities']) as f:
    viral_tags = sorted({tag_variant
                         for gene_tags in yaml.safe_load(f).values()
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

# parse experiment information
expts = pymodules.experiments.Experiments(config['experiments'],
                                          viral_tags,
                                          barcoded_viral_genes)


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
        expand(join(config['viral_tags_bcs_in_cells_dir'],
                    "{expt}_assign_viral_tags_to_cells.svg"),
               expt=expts.experiments),
        expand(join(config['viral_progeny_dir'],
                    "{expt}_viral_bc_in_progeny.csv.gz"),
               expt=expts.experiments),


# Set up report  -------------------------------------------------------------

report: 'report/workflow.rst'


# Load rules -----------------------------------------------------------------

include: 'rules/viral_tags_bcs_in_cells.smk'
include: 'rules/viral_fastq10x.smk'
include: 'rules/align_fastq10x.smk'
include: 'rules/star_refgenome.smk'
include: 'rules/fastq10x.smk'
include: 'rules/utils.smk'
include: 'rules/viral_progeny.smk'
