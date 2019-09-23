# flu sequences
This subdirectory has the flu sequences used in the experiments.

### Input Data
The [./plasmid_maps](plasmid_maps) subdirectory has the Genbank files for the Bloom lab plasmid maps for the plasmids used to generate the viruses.

#### Barcoded Segments
The plasmid maps for the barcoded segments (with `'bc'` in their Genbank file names) have some differences from the rest of the plasmids. First, these segments are carried on a pHH, rather pHW backbone. Second, these segments have a duplicated packaging signal at the 5' and 3' end of the coding sequence. Third, these segments have a 16N viral barcode inserted between the end of the coding sequence and the duplicated packaging signal.

Briefly, these plasmids were produced by amplifying the coding sequence from plasmids `2232 (HA_G155E)`, `2303 (HA_G155E-dblSyn)`, `368 (NA)`, or `2304 (NA-dblSyn)` with a primer that appends a 16N barcode to the 3' end of the CDS. They were then cloned into a backbone with an appropriate duplicated packaging signal: `2227 (HA)` or `2228 (NA)`. For complete details on generating the barcoded plasmid libraries for these segments, see David's notebook entries: [Cloning Barcoded CA09 HA](https://benchling.com/s/etr-81Sz1BB1P4xrGAhTbIF8) and [Cloning Barcoded CA09 NA](https://benchling.com/s/etr-DkKpVWAULtCRxirzKSln).

### Processing
The Python script [make_seqs_and_gtf.py](make_seqs_and_gtf.py) takes these plasmid maps and constructs files with the sequences and annotations for the rest of the analysis.
