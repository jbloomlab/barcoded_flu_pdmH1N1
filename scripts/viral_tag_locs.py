"""Locations of viral tags in 1-based indexing."""


import sys

import Bio.SeqIO

import pandas as pd


print(f"Writing output and errors to {snakemake.log.log}")
f = open(snakemake.log.log, 'w')
sys.stdout = f
sys.stderr = f
   
print(f"Parsing viral tags from {snakemake.input.viral_genbank}")
viral_tag_tups = []
for s in Bio.SeqIO.parse(snakemake.input.viral_genbank, 'genbank'):
    for f in s.features:
        if 'tag' in f.type:
            viral_tag_tups.append((s.id,
                                   f.type,
                                   int(f.location.start) + 1,
                                   int(f.location.end)))
print(f"Writing tag locations to {snakemake.output.viral_tag_locs}")
pd.DataFrame.from_records(viral_tag_tups,
                          columns=['gene', 'tag_name', 'start', 'end']
                          ).to_csv(snakemake.output.viral_tag_locs,
                                   index=False)
