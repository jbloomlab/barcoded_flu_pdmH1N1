"""Locations of viral barcodes in 1-based indexing."""


import sys

import Bio.SeqIO

import pandas as pd


print(f"Writing output and errors to {snakemake.log.log}")
f = open(snakemake.log.log, 'w')
sys.stdout = f
sys.stderr = f

print(f"Parsing viral barcodes from {snakemake.input.viral_genbank}")
viral_bc_tups = []
for s in Bio.SeqIO.parse(snakemake.input.viral_genbank, 'genbank'):
    for f in s.features:
        if f.type == 'viral_barcode':
            viral_bc_tups.append((s.id,
                                  int(f.location.start) + 1,
                                  int(f.location.end)))
print(f"Writing barcode locations to {snakemake.output.viral_bc_locs}")
pd.DataFrame.from_records(viral_bc_tups,
                          columns=['gene', 'start', 'end']
                          ).to_csv(snakemake.output.viral_bc_locs,
                                   index=False)
