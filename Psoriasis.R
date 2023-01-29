library("Seurat")
library("GEOquery")

"""gse <- getGEO("GSE151177") 
k = exprs(gse[[1]])

files <- getGEOSuppFiles("GSE151177")
dim(pData(gse[[1]]))
head(pData(gse[[1]])[, 1:3])

samples <-  gse$GSE151177_series_matrix.txt.gz@phenoData@data



gsm <- getGEO(filename=system.file("extdata/GSM4567878",package="GEOquery"))


expressionSet <- getGEO("GSE138651")
expressionSet <- expressionSet[[1]]
#You can then look at the info of the samples:
  
pData(expressionSet)
#And the assay data itself

exprs(expressionSet)
"""

seurat_obj <- ReadMtx(mtx = './GSE151177/GSM4567877_Control01_matrix.mtx.gz',
        features="./GSE151177/GSM4567877_Control01_features.tsv.gz",
        cells="./GSE151177/GSM4567877_Control01_barcodes.tsv.gz")

#cts <-  seurat_obj
seurat_obj.obj <- CreateSeuratObject(counts = seurat_obj, project = "psoriasis", min.cells = 3, min.features = 200)

seurat_obj.obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj.obj, pattern = "^MT-")
FeatureScatter(seurat_obj.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 



#VlnPlot(seurat_obj.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


seurat_obj.obj <- subset(seurat_obj.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)

# 3. Normalize data ----------
#seurat_obj.obj <- NormalizeData(seurat_obj.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
seurat_obj.obj <- NormalizeData(seurat_obj.obj)



seurat_obj.obj <- FindVariableFeatures(seurat_obj.obj, selection.method = "vst", nfeatures = 4000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj.obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_obj.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)



all.genes <- rownames(seurat_obj.obj)
seurat_obj.obj <- ScaleData(seurat_obj.obj, features = all.genes)

str(seurat_obj.obj)

# 6. Perform Linear dimensionality reduction --------------
seurat_obj.obj <- RunPCA(seurat_obj.obj, features = VariableFeatures(object = seurat_obj.obj))

# visualize PCA results
print(seurat_obj.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(seurat_obj.obj, dims = 1, cells = 500, balanced = TRUE)


# determine dimensionality of the data
ElbowPlot(seurat_obj.obj)


# 7. Clustering ------------
seurat_obj.obj <- FindNeighbors(seurat_obj.obj, dims = 1:15)

# understanding resolution
seurat_obj.obj <- FindClusters(seurat_obj.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(seurat_obj.obj@meta.data)

DimPlot(seurat_obj.obj, group.by = "RNA_snn_res.0.1", label = TRUE)


seurat_obj.obj <- RunUMAP(seurat_obj.obj, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(seurat_obj.obj, reduction = "umap")


markers <-FindAllMarkers(seurat_obj.obj,
               logfc.threshold = 0.25,
               min.pct = 0.1,
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')



