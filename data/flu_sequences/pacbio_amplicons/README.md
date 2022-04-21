# Pacbio amplicon sequences
Files in this subdirectory:

 - [fluCA09.gb](fluCA09.gb) file contains amplicon sequences used for pacbio sequencing.
 - [fluCA09_features.yaml](fluCA09_features.yaml) file contains information on how to parse the features in [fluCA09.gb](fluCA09.gb) file using [alignparse](https://jbloomlab.github.io/alignparse/).
 - [linearization_primers.tsv](linearization_primers.tsv): contains the primers used to linearize the circularized products.
 - [amplicon_to_reference.csv](amplicon_to_reference.csv): contains a lookup table to convert amplicon positions to positions in the reference CA09 sequence.

To generate amplicon templates, cDNA from the first round of PCR amplification from  chromium v3.1 single-cell sequencing library prep was used. The viral reads were first amplified using semi-specific PCR with primers aligning to a common-to-all-reads Read 1 sequence (acuired during reverse transcription in GEMs) and a unique viral segment-specific 5' termini primer that was flanked with reverse complement of Read 1. The products of semi-specific PCR  were then circularized via Read 1 site using NEB HiFi assembly kit. The cirularized DNA was then linearised and amplified using either primers aligning close to the 5' of the viral mRNA or primers aligning close to the middle of each segment (the latter approach allows to enrich for full-length segments). PCR products were then purified using Ampure beads and pooled together. Sequencing was done on PacBio SMRTcells.

To generate the GenBank amplicon sequences [flu-CA09.fasta](../flu-CA09.fasta) reference was used. [SnapGene](https://www.snapgene.com/) sequence viewer was used to stimulate the PCR procedures described above and the resulting PCR products were exported as GenBank files. Redundant sequence annotation from the GenBank file was removed manually.

The resulting GenBank files exported by SnapGene cannot be loaded directly by BioPython because of some incompatibilities in how features are annotated. To clean up and simplify the GenBank annotations, the files were loaded into [Benchling](https://www.benchling.com/) https and exported as a single GenBank file containing all amplicon targets.

[amplicon_to_reference.csv](amplicon_to_reference.csv) file contains the following columns:
- `target`: amplicon name
- `gene`: CA09 segment from which amplicon was generated
- `wt_nt`: nucleotide in CA09 at a given sited numbered in `ORF_position`
- `ORF_position`: nucleotide position in CA09
- each feature for which we are parsing mutations we want to process

[create_amplicon_to_reference.ipynb](create_amplicon_to_reference.ipynb) is the Jupyter notebook used to create [amplicon_to_reference.csv](amplicon_to_reference.csv) from the annotations in [fluCA09.gb](fluCA09.gb).

In the [fluCA09_features.yaml](fluCA09_features.yaml) NEP is defined to have only one variant tag because after splicing NEP mRNA looses variant tag, all other segments have 2 tags.

M segment primers probably don't capture M2 very well because part of the primers used to linearize M segment reads overlap the 5' splice site in the M2 and so don't align fully. M2 is, therefore, not included in [fluCA09.gb](fluCA09.gb). 
