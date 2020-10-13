Single-cell transcriptomics of cells infected with influenza virions carrying barcodes.
This experiment allows accurate detection of the number of unique virions infecting each cell and its resulting impact on the transcriptome.
The single-cell transcriptomics were performed using `10x Chromium <https://www.10xgenomics.com/solutions/single-cell/>`_.

The basic steps in the analysis are as follows:

 - FASTQ files from Illumina sequencing the 10x transcriptomic libraries are generated using `cellranger mkfastq <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq>`_.
   See `10x transcriptomics FASTQ files`_.

 - Transcriptome quantification is performed by aligning the 10x FASTQ files with STARsolo_.
   STARsolo_ also identifies the cells and creates the cell-gene matrix.

 - The viral tags and barcodes are analyzed in the aligned 10x FASTQ files.

 - The cell-gene matrix is analyzed using `scanpy <https://scanpy.readthedocs.io/>`_ and `AnnData <https://anndata.readthedocs.io/>`_.

.. _STARsolo: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
