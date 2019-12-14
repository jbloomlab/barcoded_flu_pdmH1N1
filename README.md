# Barcoded pdmH1N1 influenza virus single-cell sequencing
Single-cell sequencing of barcoded pdmH1N1 influenza virus; David Bacsik and Jesse Bloom.

## Summary of workflow and results
For a summary, see [report.html].

## Organization of repository
This repository is organized as follows.
It is based loosely on [this example snakemake repository](https://github.com/koesterlab/single-cell-rna-seq).

 - [Snakefile] is the [snakemake] file that runs the analysis.

 - [./rules/](rules) contains [snakemake] rules.

 - [config.yaml](config.yaml) contains the configuration for the analysis in [Snakefile], as described [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).

 - [cluster.yaml](cluster.yaml) contains the cluster configuration for running [Snakefile] on the Fred Hutch cluster, as described [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).

 - [./data/](data) contains the input data, specifically:

   * [./data/flu_sequences/](data/flu_sequences) gives the flu sequences used in the experiment. See the [README in that subdirectory](data/flu_sequences/README.md) for details.

   * [./data/illumina_runs_10x.csv](data/illumina_runs_10x.csv) specifies the Illumina sequencing runs of the 10X transcriptome libraries.

 - [./results/](results) is a created directory with all results, most of which are not tracked in this repository.

 - [report.html] is the report created by running [Snakefile].

## Running the analysis
To run the analysis, simple run [Snakefile] with the command:

    snakemake

And then generate [report.html] with:

    snakemake --report report.html

To do this on the Hutch cluster using [sbatch](sbatch) and the cluster configuration in [cluster.yaml](cluster.yaml), run the bash script [run_Hutch_cluster.bash](run_Hutch_cluster.bash).
You probably want to submit the script itself via [sbatch](sbatch), using:

    sbatch run_Hutch_cluster.sbatch

[report.html]: report.html
[Snakefile]: Snakefile
[snakemake]: https://snakemake.readthedocs.io
