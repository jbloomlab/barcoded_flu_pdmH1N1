#!/bin/bash

printf "Starting the linting...\n"

printf "Using 'snakemake --lint' to check Snakefile...\n"
snakemake --lint

printf "\nDoing a snakemake dry run...\n"
snakemake -n

printf "\nRunning flake8...\n"
flake8

printf "\nRunning flake8_nb...\n"
flake8_nb notebooks/*.ipynb
flake8_nb *.ipynb --ignore=E902

printf "Linting complete\n"
