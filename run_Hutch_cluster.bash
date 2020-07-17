#!/bin/bash

set -e

echo "Running snakemake..."
snakemake \
    -j 999 \
    --cluster-config cluster.yaml \
    --cluster "sbatch -c {cluster.cpus} -t {cluster.time} -J {cluster.name}" \
    --latency-wait 30 \
    --use-envmodules \
    --use-conda \
    --conda-prefix ./env \
    -R `snakemake --list-input-changes`  # https://snakemake.readthedocs.io/en/stable/project_info/faq.html#snakemake-does-not-trigger-re-runs-if-i-add-additional-input-files-what-can-i-do
echo "Run of snakemake complete."

echo "Generating report..."
snakemake --report report.html
echo "Report created."
