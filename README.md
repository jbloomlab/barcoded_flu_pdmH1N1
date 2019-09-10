# virus_hashing
Barcoded pdmH1N1 virus hashing experiment

## Experimental Outline
In this experiment, cells were infected at moderate MOI (~1) and collected at 10 hpi. I expect some cells to be uninfected, some to be infected a single virion, and some to be infected by multiple virions.

## Molecule Design
Each molecule is tagged with a 10X sample index near the Read 2 end. These indexes have four different sequences per sample.

Because the virus already had a TruSeq Read 1 primer embedded in its genome, we attempted to append a Nextera Read 1 sequence onto the end of the molecule. This should have worked every time for non-viral (host) transcripts. Those will look like this:

`P5 Adapter - Nextera Read 1 - TruSeq Read 1 - Cell Barcode - UMI - PolyA - CDS - Read 2 - 10X Sample Indices - P7`

For viral transcripts, the primer may have bound at the end of of the molecule, retaining the Cell Barcode and UMI. (For the record, this is what we want). Those molecules will look like this:

`P5 Adapter - Nextera Read 1 - TruSeq Read 1 - Cell Barcode - UMI - PolyA - TruSeq Read 1 - CDS - Read 2 - 10X Sample Indices - P7`

However, it is possible that the primer bound further into the molecule, at the second TruSeq Read 1 site, losing the Cell Barcode and UMi. Those molecules will look like this:

`P5 Adapter - Nextera Read 1 - TruSeq Read 1 - CDS - Read 2 - 10X Sample Indices - P7`


## Demultiplexing Strategy

To demultiplex these samples, I will use the 10X cellranger mkfastq software, as the indexes are unaffected by the weirdness at the Read 1 end.

The demuxed FASTQ files are stored in a relatively deep path: out/virus_hashing/outs/fastq_path/`NAME OF THE ILLUMINA RUN`/`SAMPLE`  
Each sample gets its own folder.