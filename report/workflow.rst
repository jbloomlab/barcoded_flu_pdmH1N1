Single-cell transcriptomics of cells infected with influenza virions carrying barcodes.
This experiment allows accurate detection of the number of unique virions infecting each cell and its resulting impact on the transcriptome.
The single-cell transcriptomics were performed using `10X Chromium <https://www.10xgenomics.com/solutions/single-cell/>`_.

The basic steps in the analysis are as follows:

 - The FASTQ files from Illumina sequencing the 10X transcriptomic libraries are generated using `cellranger mkfastq <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq>`_.
   See the summary results in `10X FASTQ files`_.

 - Transcriptome quantification is performed by aligning the 10X FASTQ files with STARsolo_.
   STARsolo_ also identifies the cells and creates the cell-gene matrix.
   See the summary results in `Aligning 10X FASTQs`_.

.. _STARsolo: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
