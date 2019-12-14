#!/bin/bash

set -e

echo "Running snakemake..."
snakemake \
    -j 999 \
    --cluster-config cluster.yaml \
    --cluster "sbatch -p {cluster.partition} -c {cluster.cpus} --mem={cluster.mem} -t {cluster.time} -J {cluster.name}" \
    --latency-wait 30
echo "Run of snakemake complete."

echo "Generating report..."
snakemake --report report.html
echo "Report created."
