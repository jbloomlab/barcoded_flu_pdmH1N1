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
                    "{expt}_valid_viral_barcodes_by_cell.pdf"),
               expt=expts.experiments),
        expand(join(config['viral_fastq10x_dir'],
                    "{expt}_viral_bc_by_cell_corrected.pdf"),
               expt=expts.experiments),
        expand(join(config['viral_fastq10x_dir'],
                    "{expt}_viral_bc_by_cell.svg"),
               expt=expts.experiments),
        expand(join(config['viral_fastq10x_dir'],
               "{expt}_viral_genes_by_cell.csv.gz"),
               expt=expts.experiments),
        expand(join(config['viral_fastq10x_dir'],
                    "{expt}_viral_transcript_coverage.svg"),
               expt=expts.experiments),
        expand(join(config['viral_tags_bcs_in_cells_dir'],
                    "{expt}_assign_infection_status.svg"),
               expt=expts.experiments),
        expand(join(config['viral_tags_bcs_in_cells_dir'],
                    "{expt}_assign_viral_tags_to_cells.svg"),
               expt=expts.experiments),
        expand(join(config['viral_progeny_dir'],
                    "{expt}_viral_bc_fates.svg"),
               expt=expts.expts_with_progeny_barcodes),
        expand(join(config['viral_progeny_dir'],
                    "{expt}_viral_bc_in_progeny_corrected.pdf"),
               expt=expts.expts_with_progeny_barcodes),
        expand(join(config['viral_progeny_dir'],
                    "{expt}_viral_bc_replicates.pdf"),
               expt=expts.expts_with_progeny_barcodes),
        expand(join(config['pacbio_dir'], "{expt}_ccs_summaries.svg"),
               expt=expts.expts_with_pacbio),
        expand(join(config['pacbio_dir'], "{expt}_amplicons.svg"),
               expt=expts.expts_with_pacbio),
        expand(join(config['pacbio_dir'],
                    "{expt}_plot_strand_exchange.svg"),
               expt=expts.expts_with_pacbio),
        expand(join(config['pacbio_dir'],
                    "{expt}_consensus_UMI_mutations.csv.gz"),
               expt=expts.expts_with_pacbio)


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
include: 'rules/pacbio_analysis.smk'
