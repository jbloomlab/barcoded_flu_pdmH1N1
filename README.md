# Barcoded pdmH1N1 influenza virus single-cell sequencing
Single-cell sequencing of barcoded pdmH1N1 influenza virus; David Bacsik and Jesse Bloom.

## Organization of repository
This repository is organized as follows:

 - [Snakefile](Snakefile) is the [Snakemake](https://snakemake.readthedocs.io) file that runs the analysis.

 - [./rules/](rules) contains [Snakemake](https://snakemake.readthedocs.io) rules.

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
