# Pacbio amplicon sequences

[fluCA09.gb](fluCA09.gb) file contains amplicon sequences used for pacbio sequencing. [fluCA09_features.yaml](fluCA09_features.yaml) file contains information on how to read the features in [fluCA09.gb](fluCA09.gb) file. 

To generate amplicon templates, cDNA from the first round of PCR amplification from  chromium v3.1 single-cell sequencing library prep was used. The viral reads were first amplified using semi-specific PCR with primers aligning to a common-to-all-reads Read 1 sequence (acuired during reverse transcription in GEMs) and a unique viral segment-specific 5' termini primer that was flanked with reverse complement of Read 1. The products of semi-specific PCR  were then circularized via Read 1 site using NEB HiFi assembly kit. The cirularized DNA was then linearised and amplified using either primers aligning close to the 5' of the viral mRNA or primers aligning close to the middle of each segment (the latter approach allows to enrich for full-length segments). PCR products were then purified using Ampure beads and pooled together. Sequencing was done on PacBio SMRTcells.

To generate the GenBank amplicon sequences [flu-CA09.fasta](../flu-CA09.fasta) reference was used. [SnapGene](https://www.snapgene.com/) sequence viewer was used to stimulate the PCR procedures described above and the resulting PCR products were exported as GenBank files. Redundant sequence annotation from the GenBank file was removed manually. 