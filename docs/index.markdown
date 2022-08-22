---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: default
---
### Viral transcription and progeny production in single influenza-infected cells
This page provides an interactive version of figure 5A from XXX. This plot shows the total viral transcription and progeny produced by single influenza-infected cells. Total viral transcription was calculated as the fraction of mRNA from virus in each cell. Physical progeny production (left) was calculated as the fraction of viral barcodes in the supernatant associated with each cell among barcodes assignable to any infected cell. Infectious progeny production (right) was calculated in the same way, using RNA extracted from cells infected with progeny viral supernatants.

{% include transcription_progeny_interactive.html max-width="100px" %}
<br>

### Data available on this plot
This plot shows viral transcription on the x-axis and progeny production on the y-axis. Each point represents a single influenza-infected cell. Mousing over a point provides further information:  

`cell_barcode`: Barcode that uniquely identifies each cell  
`frac_viral_UMIs`: Fraction of mRNA in the cell from frac_viral_UMIs  
`total_UMIs`: Total number of unique mRNA molecules sequenced in the cell  
`missing_viral_genes`: List of viral genes that are *not* expressed in the cell  
`mutated_genes`: List of viral genes that have any mutation  
`viral_mutations`: List of mutations identified in the viral genes and classification of the mutation (e.g. synonymous, non-synonymous, non-coding)
`mutation_support`: Number of UMIs sequenced in the PacBio data supporting a consensus mutation  

## Source data and code
The input data and code used to generate this plot are in the [final_analysis.py.ipynb](https://github.com/jbloomlab/barcoded_flu_pdmH1N1/blob/main/final_analysis.py.ipynb) notebook.
