# Barcoded pdmH1N1 influenza virus single-cell sequencing
Single-cell sequencing of barcoded pdmH1N1 influenza virus; David Bacsik and Jesse Bloom.

## Summary of workflow and results
For a summary, see the `report.html` file that is placed in the `./results/` subdirectory.

## Organization of repository
This repository is organized as followed (based loosely on [this example snakemake repository](https://github.com/koesterlab/single-cell-rna-seq)):

 - [Snakefile] is the [snakemake] file that runs the analysis.

 - [environment.yml](environment.yml) and [environment_unpinned.yml](environment_unpinned.yml) give the version pinned and unpinned [conda] environment for the analysis.

 - [config.yaml](config.yaml) contains the configuration for the analysis.

 - [cluster.yaml](cluster.yaml) contains the cluster configuration for running tha analysis on the Fred Hutch cluster.

 - [./rules/](rules) contains [snakemake] rules.

 - [./notebooks/](notebooks) contains [Jupyter notebooks](https://jupyter.org/) that are run by [Snakefile] using the [snakemake notebook functionality](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#jupyter-notebook-integration).

 - [./scripts/](scripts) contains scripts used by [Snakefile].

 - [./pymodules/](pymodules) contains Python modules with some functions used by [Snakefile].

 - [./report/](report) contains workflow description and captions used to create the [snakemake report].

 - [./data/](data) contains the input data, specifically:

   * [./data/flu_sequences/](data/flu_sequences) gives the flu sequences used in the experiment. See the [README in that subdirectory](data/flu_sequences/README.md) for details.

   * [./data/flu_sequences/pacbio_amplions](data/flu_sequences/pacbio_amplions) gives the famplicon sequences generated for pacbio sequencing. See the [README in that subdirectory](data/flu_sequences/pacbio_amplicons/README.md) for details.

 - [./results/](results) is a created directory with all results, most of which are not tracked in this repository.


## Running the analysis

### Installing software
The [conda] environment for this repo is specified in [environment.yml](environment.yml); note also that an **unpinned** version of this environment is specified in [environment_unpinned.yml](environment_unpinned.yml).
If you are on the Hutch cluster and set up to use the *BloomLab* [conda] installation, then this environment is already built and you can activate it simply with:

    conda activate barcoded_flu_pdmH1N1

Otherwise you need to first build the [conda] environment from [environment.yml](environment.yml) and then activate it as above.

In addition to building and activating the [conda] environment, you also need to install [cellranger] and [bcl2fastq] into the current path; the current analysis uses [cellranger] version 4.0.0 and [bcl2fastq] version 2.20.

### Run the analysis
Once the *barcoded_flu_pdmH1N1* [conda] environment and other software have been activated, simply enter the commands to run [Snakefile] and then generate a [snakemake report], at `./results/report.html`.
These commands with the configuration for the Fred Hutch cluster are in the shell script. [run_Hutch_cluster.bash](run_Hutch_cluster.bash).
You probably want to submit the script itself via [sbatch](sbatch), using:

    sbatch run_Hutch_cluster.sbatch

When the Snakeamke pipeline has run completely, the processed output data is exported to a CSV file at `results/viral_fastq10x/{expt}_integrate_data.csv`. This CSV file is used to perform the final analysis and generate figures in the `final_analysis.py.ipynb` notebook.

## Linting the code
Before you commit a new branch, you should run the linting in [lint.bash](lint.bash) with the command:

    bash ./lint.bash

This script runs:

 - [snakemake linting](https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html#best-practices)

 - a [snakemake dry run](https://snakemake.readthedocs.io/en/stable/project_info/faq.html#my-workflow-is-very-large-how-do-i-stop-snakemake-from-printing-all-this-rule-job-information-in-a-dry-run)

 - a [flake8](https://flake8.pycqa.org/) analysis of the Python code

 - a [flake8_nb](https://flake8-nb.readthedocs.io/) analysis of the Jupyter notebooks.

For the Jupyter notebook linting, it may be easiest to lint while you are still developing notebook with run cells rather then before you put the empty notebook in [./notebooks/](notebooks), as the linting results are labeled by cell run number.

[Snakefile]: Snakefile
[snakemake]: https://snakemake.readthedocs.io
[snakemake report]: https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html
[conda]: https://docs.conda.io/projects/conda/en/latest/index.html
[cellranger]: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
[bcl2fastq]: https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html
