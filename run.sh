#!/bin/bash
#
#SBATCH --job-name=demux
#SBATCH --output=demux_results.txt
#
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=32000

srun python make_samplefile_forcellranger.py
srun cellranger mkfastq --id=virus_hashing \
                     --run=/shared/ngs/illumina/bloom_lab/190830_M03100_0474_000000000-CL3WR \
                     --csv=simple-samples.csv