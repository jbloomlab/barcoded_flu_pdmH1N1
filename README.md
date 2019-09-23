# Barcoded pdmH1N1 influenza virus single-cell sequencing
Single-cell sequencing of barcoded pdmH1N1 influenza virus; David Bacsik and Jesse Bloom.

## Organization of repository
This repository is organized as follows:

 - [Snakefile](Snakefile) is the [Snakemake](https://snakemake.readthedocs.io) file that runs the analysis.

 - [./pymodules/](pymodules) contains custom Python modules used in analysis.

 - [config.yaml](config.yaml) contains the configuration for the analysis in [Snakefile](Snakefile), as described [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).

 - [cluster.yaml](cluster.yaml) contains the cluster configuration for running [Snakefile](Snakefile) on the Fred Hutch cluster, as described [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).

 - [./data/](data) contains the input data, specifically:

   * [./data/flu_sequences/](data/flu_sequences) gives the flu sequences used in the experiment. See the [README in that subdirectory](data/flu_sequences/README.md) for details.

   * [./data/illumina_runs_10x.csv](data/illumina_runs_10x.csv) specifies the Illumina sequencing runs of the 10X transcriptome libraries.

 - [./results/](results) is a created directory with all results, most of which are not tracked in this repository.

## Running the analysis
To run the analysis, simple run [Snakefile](Snakefile) with the command:

    snakemake

To do this on the Hutch cluster using [sbatch](sbatch) and the cluster configuration in [cluster.yaml](cluster.yaml), run the bash script [run_Huch_cluster.bash](run_Hutch_cluster.bash).
You probably want to submit the script itself via [sbatch](sbatch), using:

    sbatch run_Hutch_cluster.sbatch

## Samples analyzed
Here are brief descriptions of the samples analyzed here.

### dual_truseq_pilot
MDCK cells were infected with pdmH1N1 flu at MOI of about 1 and collected 10 hpi.
The virus contains a barcode followed by a TruSeq read 1 primer embedded in its genome in HA and NA segments.
The cell transcriptomes were analyzed using the 10X Chromium workflow (v2 reagents).
But since the virus already has a TruSeq read 1 primer site, the standard 10X SI primer was replaced with a custom primer that appended the Nextera read 1 primer downstream of the TruSeq binding site, and the libraries were sequenced using the Nextera primer.
Therefore, the desired sequence reading in reverse from the 3' end of the tagged cDNA molecules is:

    P5 Adapter - Nextera Read 1 - TruSeq Read 1 - Cell Barcode - UMI - PolyA - CDS -  Read 2 - 10X Sample Indices - P7

However, since the viral HA and NA genes also have a TruSeq primer binding site, the custom primer could also bind there, yielding molecules that look like this: 

    P5 Adapter - Nextera Read 1 - TruSeq Read 1 - Viral Barcode - Viral CDS - Read 2 - 10X Sample Indices - P7

The custom primer was added at 100 uM and 10 uM, corresponding to *pdmH1N1_dual_truseq_pilot_100uM* and *pdmH1N1_dual_truseq_pilot_10uM*.
