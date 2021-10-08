#!/usr/bin/env Rscript

# Expected arguments
# args[1]: filtered MTX file
# args[2]: filtered cell barcode file
# args[3]: filtered feature name file
# args[4]: raw MTX file
# args[5]: raw cell barcode file
# args[6]: raw feature name file
# args[7]: corrected MTX file

args <- commandArgs(trailingOnly = TRUE)

library(SoupX)
library(Seurat)
library(Matrix)

toc <- ReadMtx(
  mtx = args[1], 
  cells = args[2],
  features = args[3]
)
tod <- ReadMtx(
  mtx = args[4], 
  cells = args[5],
  features = args[6]
)

# setup SoupX object
sc <- SoupChannel(tod = tod, toc = toc)

# use flu genes to estimate contamination fraction
flu.genes <- grep(pattern = "flu", x = rownames(sc$soupProfile), value = TRUE)
useToEst <- estimateNonExpressingCells(sc, nonExpressedGeneList = list(flu = flu.genes), clusters = FALSE)
sc <- calculateContaminationFraction(sc, list(flu = flu.genes), useToEst = useToEst)


# Rough clustering of cells (speeds up correction)
ob <- CreateSeuratObject(counts = toc)
ob <- NormalizeData(object = ob)
ob <- FindVariableFeatures(object = ob)
ob <- ScaleData(object = ob)
ob <- RunPCA(object = ob)
ob <- FindNeighbors(object = ob, dims = 1:10)
ob <- FindClusters(object = ob, resolution = 0.3)
sc <- setClusters(sc = sc, clusters = setNames(object = ob$seurat_clusters, nm = Cells(x = ob)))

# Correct expression
out <- adjustCounts(sc = sc)
writeMM(obj = out, file = args[7])
