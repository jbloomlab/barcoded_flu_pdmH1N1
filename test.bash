#!/bin/bash

printf "Using 'snakemake --lint' to check Snakefile...\n"
snakemake --lint

printf "\nDoing a snakemake dry run...\n"
snakemake -n

printf "\nRunning flake8...\n"
flake8
