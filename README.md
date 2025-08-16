# scRNA-seq Analysis Tutorial (R/Seurat, step-by-step)
### *For the Chinese version of this guide, please see the **[中文版本](./Tutorials/README_CN.md)**.*
---

Author: Charlene  
Version: v1.0  
Last updated: 2025-08-16

> This tutorial targets beginners and researchers who need a fast, practical workflow centered on **R + Seurat**. It covers installation and basic syntax, the Seurat object structure, QC, normalization, dimensionality reduction & clustering, annotation, doublet detection, integration, DE/GSEA, subsetting & reclustering, and reproducibility. A small example using `pbmc_small` (lightweight, no download required) is included, plus advanced sections that require external data (enable as needed).

---

## Table of Contents
- [0. Setup & Environment](#0-setup--environment)
- [1. Basic Syntax & Getting Help](#1-basic-syntax--getting-help)
- [2. Seurat Object & Data Structure (Core)](#2-seurat-object--data-structure-core)
- [3. Quick Start (pbmc_small)](#3-quick-start-pbmc_small)
- [4. Loading Data (10x/Matrix/RDS)](#4-loading-data-10xmatrixrds)
- [5. Quality Control (QC) & Filtering](#5-quality-control-qc--filtering)
- [6. Normalization & Feature Selection: SCTransform vs LogNormalize](#6-normalization--feature-selection-sctransform-vs-lognormalize)
- [7. Dimensionality Reduction, Neighbor Graph, Clustering & UMAP](#7-dimensionality-reduction-neighbor-graph-clustering--umap)
- [8. Markers & Initial Manual Annotation](#8-markers--initial-manual-annotation)
- [9. Automatic Annotation (SingleR, optional)](#9-automatic-annotation-singler-optional)
- [10. Doublet Detection (scDblFinder)](#10-doublet-detection-scdblfinder)
- [11. Multi-sample / Batch Integration (SCT / RPCA)](#11-multi-sample--batch-integration-sct--rpca)
- [12. Cell Cycle & Regression](#12-cell-cycle--regression)
- [13. Differential Expression (DE) & GSEA](#13-differential-expression-de--gsea)
- [14. Subsetting & Reclustering](#14-subsetting--reclustering)
- [15. Export & Interoperability](#15-export--interoperability)
- [16. Reproducibility & Common Pitfalls](#16-reproducibility--common-pitfalls)
- [Appendix A: Minimal Pipeline Template](#appendix-a-minimal-pipeline-template)
- [Appendix B: Glossary](#appendix-b-glossary)

---

## 0. Setup & Environment

- Recommended: R ≥ 4.3 and RStudio Desktop.  
- Text encoding: Use UTF‑8 end-to-end. On Windows, you can add `locale = locale(encoding = "UTF-8")` (from `readr`) when reading/writing files.

**Install/load common packages (install once; use `library()` at runtime):**
```r
install.packages(c(
  "Seurat","SeuratData","patchwork","tidyverse","remotes",
  "rmarkdown","knitr","Matrix"
))

# Bioconductor packages (as needed)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "SingleR","celldex","scDblFinder","DropletUtils",
  "clusterProfiler","fgsea","org.Hs.eg.db","org.Mm.eg.db"
))
install.packages(c("msigdbr","SeuratDisk"))

# Optional: speed up SCTransform
BiocManager::install("glmGamPoi")

# Load when needed
library(Seurat); library(SeuratData); library(patchwork); library(tidyverse)
```

---

## 1. Basic Syntax & Getting Help

**Assign & print**
```r
x <- 3 * 7; x  # 21
```

**Help**
```r
?NormalizeData            # or help("NormalizeData")
help.search("integration")
apropos("FindCluster")    # fuzzy search
```

**Pipes (base R 4.1+ `|>` )**
```r
1:5 |> mean()
mtcars |> head(3) |> summary()
```

**Tips:** Use `class()`, `str()`, `dim()`, `head()`, `summary()` to inspect objects; after an error, check `traceback()`.

---

## 2. Seurat Object & Data Structure (Core)

### 2.1 Skeleton
- **Assays** (`obj@assays`): data layers (RNA / SCT / integrated / ADT …).  
- **meta.data** (`obj@meta.data`): per-cell metadata (rows = cells).  
- **reductions** (`obj@reductions`): embeddings (`pca`, `umap`, `tsne`, `lsi`).  
- **graphs** (`obj@graphs`): KNN/SNN/WNN graphs.  
- **commands** (`obj@commands`): pipeline steps & parameters.  
- **tools/misc**: extra analysis or scratch space.  

**Best practice:** Prefer accessors (`Assays()`, `DefaultAssay()`, `GetAssayData()`, `Idents()`, `Embeddings()`, `Cells()`, `Features()`) over `@` whenever possible.

### 2.2 Assay slots & variable features
- `counts` (raw counts, sparse, integers)  
- `data` (normalized expression, e.g., `log1p(CPM)`)  
- `scale.data` (scaled expression used for PCA/plots)  
- `VariableFeatures(obj)` (variable genes)

```r
obj <- pbmc_small
DefaultAssay(obj) <- "RNA"

m_counts <- GetAssayData(obj, slot = "counts")
m_data   <- GetAssayData(obj, slot = "data")
dim(m_counts); dim(m_data)
head(VariableFeatures(obj))
```

### 2.3 `meta.data` and `Idents`
```r
head(obj@meta.data, 2)
obj$dummy_tag <- sample(c("A","B"), ncol(obj), TRUE)  # quickly add a column
table(obj$dummy_tag)

Idents(obj) <- "seurat_clusters"  # set active identity
```

### 2.4 DimReduc: PCA/UMAP internals
- `@cell.embeddings`: coordinates (cells × dimensions)  
- `@feature.loadings`: feature loadings (method-dependent)  
- `@stdev`: PCA standard deviations  
- `@key`: coordinate prefix (e.g., `PC_` / `UMAP_`)  

```r
obj <- NormalizeData(obj) |>
  FindVariableFeatures(nfeatures = 200) |>
  ScaleData()

obj <- RunPCA(obj, verbose = FALSE) |>
  RunUMAP(dims = 1:10)

head(Embeddings(obj, "pca")[,1:3])
```

### 2.5 Graphs: KNN/SNN/WNN
```r
obj <- FindNeighbors(obj, dims = 1:10)
names(obj@graphs)  # e.g., "RNA_snn"
```

### 2.6 commands / tools / misc
```r
names(obj@commands)[1:5]           # see major steps & parameters
str(obj@tools, max.level = 1)
str(obj@misc,  max.level = 1)
```

### 2.7 Cheat Sheet: “Where is X stored?”
| Need               | Where                | How to access                                      |
|--------------------|----------------------|----------------------------------------------------|
| Raw counts         | `Assay@counts`       | `GetAssayData(obj, slot = "counts")`               |
| Normalized expr    | `Assay@data`         | `GetAssayData(obj, slot = "data")`                 |
| Scaled expr        | `Assay@scale.data`   | `GetAssayData(obj, slot = "scale.data")`           |
| Variable genes     | Assay/object         | `VariableFeatures(obj)`                            |
| Cell metadata      | `obj@meta.data`      | `obj$column` or `obj@meta.data`                    |
| Clusters/identity  | Idents               | `Idents(obj)` / `RenameIdents()`                   |
| PCA/UMAP coords    | DimReduc             | `Embeddings(obj, "pca"/"umap")`                    |
| KNN/SNN graphs     | `obj@graphs`         | after `FindNeighbors()`, check `names(obj@graphs)` |
| Pipeline record    | `obj@commands`       | `names(obj@commands)`                              |
| Cell/Gene names    | Cells/Features       | `Cells(obj)` / `Features(obj)`                     |

---

## 3. Quick Start (pbmc_small)

Goal: run the minimal pipeline end-to-end; build intuition for the Seurat object.

```r
library(Seurat); library(patchwork)

obj <- pbmc_small
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 200)
obj <- ScaleData(obj)
obj <- RunPCA(obj, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:10)

DimPlot(obj, reduction = "umap", label = TRUE) + NoLegend()
FeaturePlot(obj, features = head(VariableFeatures(obj), 4))
```

---

## 4. Loading Data (10x/Matrix/RDS)

```r
# 10x folder (contains barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz)
# counts <- Read10X(data.dir = "data/10x_pbmc/")
# obj <- CreateSeuratObject(counts, project = "PBMC", min.cells = 3, min.features = 200)

# RDS (recommended for saving/loading Seurat objects)
# saveRDS(obj, "data/seurat_obj.rds")
# obj <- readRDS("data/seurat_obj.rds")
```

---

## 5. Quality Control (QC) & Filtering

```r
mt_pat <- if (any(grepl("^MT-", rownames(obj)))) "^MT-" else "^mt-"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pat)

VlnPlot(obj, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(obj, feature1 = "nCount_RNA",   feature2 = "nFeature_RNA") +
FeatureScatter(obj, feature1 = "nCount_RNA",   feature2 = "percent.mt")

# Use MAD as a rough guide (adjust to your data)
mad_cut <- function(x, nmads=3){
  med <- median(x); m <- mad(x)
  c(max(0, med - nmads*m), med + nmads*m)
}
# lim_feat  <- mad_cut(obj$nFeature_RNA)
# lim_count <- mad_cut(obj$nCount_RNA)
# mt_upper  <- 20  # Humans often < 10–20%
# obj <- subset(obj, subset =
#   nFeature_RNA > lim_feat[1] & nFeature_RNA < lim_feat[2] &
#   nCount_RNA   > lim_count[1] & nCount_RNA   < lim_count[2] &
#   percent.mt   < mt_upper)
```

> Heuristics (tune by dataset): low `nFeature_RNA` suggests empty droplets; high `nFeature_RNA` suggests doublets; for humans, `percent.mt` is often < 10–20%.  
> Advanced: remove empty droplets with `DropletUtils::emptyDrops`; correct ambient RNA with `SoupX`.

---

## 6. Normalization & Feature Selection: SCTransform vs LogNormalize

**Recommendation: `SCTransform()`** (stable; can regress `percent.mt` / cell-cycle scores).
```r
# SCTransform (recommended)
# obj <- SCTransform(obj, variable.features.n = 3000, vars.to.regress = "percent.mt", verbose = FALSE)
# DefaultAssay(obj) <- "SCT"
```

**Traditional workflow (good for teaching/small data):**
```r
# LogNormalize
# obj <- NormalizeData(obj)
# obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
# obj <- ScaleData(obj, vars.to.regress = "percent.mt")
# DefaultAssay(obj) <- "RNA"
```

---

## 7. Dimensionality Reduction, Neighbor Graph, Clustering & UMAP

```r
# obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
# ElbowPlot(obj, ndims = 50)
# dims_to_use <- 1:30
# obj <- FindNeighbors(obj, dims = dims_to_use)
# obj <- FindClusters(obj, resolution = 0.4)
# obj <- RunUMAP(obj, dims = dims_to_use)
# DimPlot(obj, reduction = "umap", label = TRUE) + NoLegend()
```

---

## 8. Markers & Initial Manual Annotation

```r
# DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else "RNA"
# Idents(obj) <- "seurat_clusters"
# markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# head(markers)
# VlnPlot(obj, features = c("MS4A1","LYZ"), pt.size = 0)
# FeaturePlot(obj, features = c("MS4A1","CD3D","LYZ","GNLY"))
```

---

## 9. Automatic Annotation (SingleR, optional)

```r
# library(SingleR); library(celldex); library(SingleCellExperiment)
# ref <- celldex::HumanPrimaryCellAtlasData()
# sce <- as.SingleCellExperiment(obj)
# pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
# obj$SingleR_label <- pred$labels
# DimPlot(obj, group.by = "SingleR_label", label = TRUE, repel = TRUE)
```

---

## 10. Doublet Detection (scDblFinder)

```r
# library(scDblFinder); library(SingleCellExperiment)
# sce <- as.SingleCellExperiment(obj)
# sce <- scDblFinder(sce)
# table(sce$scDblFinder.class)
# obj$doublet <- sce$scDblFinder.class
# obj <- subset(obj, subset = doublet == "singlet")
```

---

## 11. Multi-sample / Batch Integration (SCT / RPCA)

**SCT integration (recommended)**
```r
# features  <- SelectIntegrationFeatures(object.list = list_objs, nfeatures = 3000)
# list_objs <- PrepSCTIntegration(object.list = list_objs, anchor.features = features, verbose = FALSE)
# anchors   <- FindIntegrationAnchors(object.list = list_objs, normalization.method = "SCT",
#                                     anchor.features = features, verbose = FALSE)
# integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
# DefaultAssay(integrated) <- "integrated"
```

**RPCA integration (good for similar compositions / larger datasets)**
```r
# list_objs <- lapply(list_objs, \(x){
#   x |> NormalizeData() |> FindVariableFeatures() |> ScaleData() |> RunPCA(verbose = FALSE)
# })
# features <- SelectIntegrationFeatures(object.list = list_objs, nfeatures = 3000)
# anchors  <- FindIntegrationAnchors(object.list = list_objs, anchor.features = features,
#                                    reduction = "rpca", dims = 1:30)
# integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
```

---

## 12. Cell Cycle & Regression

```r
# s.genes   <- Seurat::cc.genes.updated.2019$s.genes
# g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
# obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# obj <- SCTransform(obj, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = FALSE)
```

---

## 13. Differential Expression (DE) & GSEA

**Differential expression**
```r
# Idents(obj) <- "seurat_clusters"
# de_0_vs_all <- FindMarkers(obj, ident.1 = 0, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
# head(de_0_vs_all)
```

**GSEA (`clusterProfiler` + `fgsea` + `msigdbr`)**
```r
# library(clusterProfiler); library(msigdbr); library(fgsea)
# species <- "Homo sapiens"
# msig <- msigdbr(species = species, category = "H") |> dplyr::select(gs_name, gene_symbol)
# path_list <- split(msig$gene_symbol, msig$gs_name)
# mk <- FindMarkers(obj, ident.1 = 0, min.pct = 0.1, logfc.threshold = 0.1)
# rk <- mk$avg_log2FC; names(rk) <- rownames(mk); rk <- sort(rk, decreasing = TRUE)
# fg <- fgsea(pathways = path_list, stats = rk, minSize = 15, maxSize = 500, nperm = 1000)
```

> Common error: “`No gene can be mapped`” → check **species**, **ID type** (`SYMBOL`/`ENTREZID`/`ENSEMBL`), and **case** (humans: upper case; mouse: lower case / capitalization).

---

## 14. Subsetting & Reclustering

```r
# tcell_genes <- c("CD3D","CD3E","CD2")
# obj$TcellFlag <- Matrix::colSums(GetAssayData(obj, slot = "data")[tcell_genes, , drop = FALSE]) > 0
# tcell <- subset(obj, subset = TcellFlag)
# tcell <- SCTransform(tcell, variable.features.n = 2000, verbose = FALSE)
# tcell <- RunPCA(tcell); tcell <- FindNeighbors(tcell, dims = 1:20)
# tcell <- FindClusters(tcell, resolution = 0.6); tcell <- RunUMAP(tcell, dims = 1:20)
```

---

## 15. Export & Interoperability

```r
# saveRDS(obj, "export/seurat_obj.rds")

# library(SeuratDisk)
# SaveH5Seurat(obj, filename = "export/seurat_obj.h5seurat", overwrite = TRUE)
# Convert("export/seurat_obj.h5seurat", dest = "h5ad", overwrite = TRUE)  # for Scanpy
```

---

## 16. Reproducibility & Common Pitfalls

```r
# if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
# renv::init(); renv::snapshot()
# sessionInfo()
```

**Common pitfalls**
- Mito/ribo prefixes: Humans `^MT-` / `^RP[LS]`; Mouse `^mt-` / `^Rp[ls]`.  
- Over/under-filtering: inspect distributions first; MAD is only a heuristic.  
- Residual batch effects after integration: tune `nfeatures`/`dims`, try Harmony/fastMNN; remove extreme clusters first.  
- Doublets: increase expected rate or retune; look for high `nCount_RNA` with multiple lineage markers.  
- GSEA mapping: ensure species/ID/case match; handle duplicate symbols.  

---

## Appendix A: Minimal Pipeline Template

```r
# library(Seurat); set.seed(123)
# counts <- Read10X("data/10x/")
# obj <- CreateSeuratObject(counts, min.cells = 3, min.features = 200)
# obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = if (any(grepl("^MT-", rownames(obj)))) "^MT-" else "^mt-")
# obj <- SCTransform(obj, variable.features.n = 3000, vars.to.regress = "percent.mt", verbose = FALSE)
# obj <- RunPCA(obj); obj <- FindNeighbors(obj, dims = 1:30)
# obj <- FindClusters(obj, resolution = 0.6); obj <- RunUMAP(obj, dims = 1:30)
# DimPlot(obj, label = TRUE)
```

---

## Appendix B: Glossary

- **Feature/Gene**: gene (rows); **Cell/Barcode**: cells (columns).  
- **Assay**: data layer (e.g., RNA, SCT, integrated).  
- **Reduction**: dimensionality reduction (PCA/UMAP/TSNE).  
- **Graph**: cell–cell neighbor relationships (KNN/SNN/WNN).  
- **Anchors**: correspondences used in integration/label transfer.  
- **Doublet**: two cells mixed under the same barcode.  

---

License: This tutorial is released under **CC BY 4.0**.
