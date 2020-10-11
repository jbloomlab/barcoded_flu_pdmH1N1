# Barcoded pdmH1N1 influenza virus single-cell sequencing
Single-cell sequencing of barcoded pdmH1N1 influenza virus; David Bacsik and Jesse Bloom.

## Summary of workflow and results
For a summary, see [report.html].

## Organization of repository
This repository is organized as followed (based loosely on [this example snakemake repository](https://github.com/koesterlab/single-cell-rna-seq)):

 - [Snakefile] is the [snakemake] file that runs the analysis.

 - [config.yaml](config.yaml) contains the configuration for the analysis in [Snakefile], as described [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).

 - [cluster.yaml](cluster.yaml) contains the cluster configuration for running [Snakefile] on the Fred Hutch cluster, as described [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).

 - [./rules/](rules) contains [snakemake] rules.

 - [./notebooks/](notebooks) contains [Jupyter notebooks](https://jupyter.org/) that are run by [Snakefile] using [papermill parameterization](https://papermill.readthedocs.io/).

 - [./pymodules/](pymodules) contains Python modules with some functions used by [Snakefile].

 - [./report/](report) contains workflow description and captions used to create the [snakemake report].

 - [./data/](data) contains the input data, specifically:

   * [./data/flu_sequences/](data/flu_sequences) gives the flu sequences used in the experiment. See the [README in that subdirectory](data/flu_sequences/README.md) for details.

   * [./data/illumina_runs_10x.csv](data/illumina_runs_10x.csv) specifies the Illumina sequencing runs of the 10X transcriptome libraries.

 - [report.html] is the report created by running [Snakefile].

 - [./results/](results) is a created directory with all results, most of which are not tracked in this repository.


## Running the analysis
The [conda] environment for this repo is specified in [environment.yml](environment.yml); note also that an **unpinned** version of this environment is  ecified[environment_unpinned.yml](environment_unpinned.yml).
If you are on the Hutch cluster and set up to use the *BloomLab* [conda] installation, then this environment is already built and you can activate it simply with:

    conda activate barcoded_flu_pdmH1N1

Otherwise you need to first build the [conda] environment from [environment.yml](environment.yml) and then activate it as above.

Once the *barcoded_flu_pdmH1N1* [conda] environment has been activated, simply run [Snakefile] with the command:

    snakemake

And then generate the [snakemake report], [report.html], with:

    snakemake --report report.html

To run on the Hutch cluster using [sbatch](sbatch) and the cluster configuration in [cluster.yaml](cluster.yaml), run the bash script [run_Hutch_cluster.bash](run_Hutch_cluster.bash).
You probably want to submit the script itself via [sbatch](sbatch), using:

    sbatch run_Hutch_cluster.sbatch

[report.html]: report.html
[Snakefile]: Snakefile
[snakemake]: https://snakemake.readthedocs.io
[snakemake report]: https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html
[conda]: https://docs.conda.io/projects/conda/en/latest/index.html


## Experimental Details
Generation of the virus libraries used in the experiments here are described [here](https://benchling.com/s/etr-tP5BHW8xCbkT3AKxWkc2). Expansion of these virus libraries is described [here](https://benchling.com/s/etr-ywy1Wbt7qcLXYnj5jhhs).

The following experiments are analyzed in this repository:
* [wt_rapidpilot](https://benchling.com/s/etr-Q28fCd1kprRNxAd0v5Hg): This experiment was performed using a small-scale rescue of the WT virus.
* [hashing_trial1](https://benchling.com/s/etr-i9I0yHiFb0P8wHCxosim): This experiment was performed using the standard WT and dblSyn virus libraries. Infection volume was chosen based on HA expression measured by flow cytometry.
* [hashing_trial2](https://benchling.com/s/etr-W8urOmOAQ7L6U4HAXMNy): This is experiment also used the standard WT and dblSyn virus libraries. Infection volume was chosen based on the results of `hashing_trial1` and flow cytometery. The innoculum volume for WT was about 12-fold higher, and the innoculum volume for dblSyn was about 24-fold higher.
