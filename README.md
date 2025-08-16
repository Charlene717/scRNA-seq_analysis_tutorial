# scRNA-seq 分析教程（R/Seurat，由淺入深）
作者：Charlene  
版本：v1.2  
最後更新：2025-08-16

> 本教程面向初學者與需要快速上手的研究者，以 **R + Seurat** 為核心，從安裝與基礎語法，到 Seurat 物件結構、QC、正規化、降維分群、註解、雙重體偵測、整合、DE/GSEA、子群重分析與可重現性，完整呈現一條可實作的分析路徑。內含可立即執行的 `pbmc_small` 範例（輕量且無需下載），以及需要外部資料的進階段落（請依需求啟用）。

---

## 目錄
- [0. 安裝與環境](#0-安裝與環境)
- [1. 基本語法與求助](#1-基本語法與求助)
- [2. Seurat 物件與資料結構（重點）](#2-seurat-物件與資料結構重點)
- [3. 快速上手（pbmc_small）](#3-快速上手pbmc_small)
- [4. 載入資料（10x/Matrix/RDS）](#4-載入資料10xmatrixrds)
- [5. 品質管制（QC）與過濾](#5-品質管制qc與過濾)
- [6. 正規化與特徵選擇：SCTransform vs LogNormalize](#6-正規化與特徵選擇sctransform-vs-lognormalize)
- [7. 降維、鄰居圖、分群與 UMAP](#7-降維鄰居圖分群與-umap)
- [8. Marker 與初步註解（手動）](#8-marker-與初步註解手動)
- [9. 自動註解（SingleR，可選）](#9-自動註解singler可選)
- [10. 雙重體偵測（scDblFinder）](#10-雙重體偵測scdblfinder)
- [11. 多樣本/批次整合（SCT / RPCA）](#11-多樣本批次整合sct--rpca)
- [12. 細胞週期與回歸](#12-細胞週期與回歸)
- [13. 差異表達（DE）與 GSEA](#13-差異表達de-與-gsea)
- [14. 亞群重分析（Subsetting & Reclustering）](#14-亞群重分析subsetting--reclustering)
- [15. 資料導出與跨工具互通](#15-資料導出與跨工具互通)
- [16. 可重現性與常見陷阱](#16-可重現性與常見陷阱)
- [附錄A：Minimal Pipeline 模板](#附錄aminimal-pipeline-模板)
- [附錄B：名詞對照（Glossary）](#附錄b名詞對照glossary)

---

## 0. 安裝與環境

- 建議環境：R ≥ 4.3 與 RStudio Desktop。  
- 文字編碼：建議全程使用 UTF-8；Windows 讀寫可加 `locale = locale(encoding = "UTF-8")`（`readr`）。

**安裝/載入常用套件（安裝一次；使用時以 `library()` 載入）：**
```r
install.packages(c(
  "Seurat","SeuratData","patchwork","tidyverse","remotes",
  "rmarkdown","knitr","Matrix"
))

# Bioconductor 套件（視分析需求）
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "SingleR","celldex","scDblFinder","DropletUtils",
  "clusterProfiler","fgsea","org.Hs.eg.db","org.Mm.eg.db"
))
install.packages(c("msigdbr","SeuratDisk"))

# 可選：SCTransform 加速
BiocManager::install("glmGamPoi")

# 需用時載入
library(Seurat); library(SeuratData); library(patchwork); library(tidyverse)
```

---

## 1. 基本語法與求助

**指派與列印**
```r
x <- 3 * 7; x  # 21
```

**求助**
```r
?NormalizeData            # 或 help("NormalizeData")
help.search("integration")
apropos("FindCluster")    # 模糊搜
```

**管線（R 4.1+ 原生 `|>`）**
```r
1:5 |> mean()
mtcars |> head(3) |> summary()
```

**小技巧**：常用檢視函數 `class()`、`str()`、`dim()`、`head()` 與 `summary()`；報錯後可看 `traceback()`。

---

## 2. Seurat 物件與資料結構（重點）

### 2.1 大骨架
- **Assays**（`obj@assays`）：資料層（RNA / SCT / integrated / ADT…）。  
- **meta.data**（`obj@meta.data`）：細胞層級中介資料（列 = 細胞）。  
- **reductions**（`obj@reductions`）：降維結果（`pca`、`umap`、`tsne`、`lsi`）。  
- **graphs**（`obj@graphs`）：KNN/SNN/WNN 圖。  
- **commands**（`obj@commands`）：流程與參數紀錄。  
- **tools/misc**：延伸分析或自訂暫存。  

**最佳實踐**：優先用存取器（`Assays()`、`DefaultAssay()`、`GetAssayData()`、`Idents()`、`Embeddings()`、`Cells()`、`Features()`）而非 `@`。

### 2.2 Assay 三大槽位與高變異基因
- `counts`（原始計數，稀疏矩陣，整數）  
- `data`（正規化表達，如 `log1p(CPM)`）  
- `scale.data`（縮放後表達，多用於 PCA/圖形）  
- `VariableFeatures(obj)`：高變異基因

```r
obj <- pbmc_small
DefaultAssay(obj) <- "RNA"

m_counts <- GetAssayData(obj, slot = "counts")
m_data   <- GetAssayData(obj, slot = "data")
dim(m_counts); dim(m_data)
head(VariableFeatures(obj))
```

### 2.3 `meta.data` 與 `Idents`
```r
head(obj@meta.data, 2)
obj$dummy_tag <- sample(c("A","B"), ncol(obj), TRUE)  # 快速新增欄位
table(obj$dummy_tag)

Idents(obj) <- "seurat_clusters"  # 設定活躍身分
```

### 2.4 DimReduc：PCA/UMAP 結構
- `@cell.embeddings`：細胞座標（細胞 × 維度）  
- `@feature.loadings`：特徵負載（方法依據）  
- `@stdev`：PCA 標準差  
- `@key`：座標前綴（如 `PC_`/`UMAP_`）  

```r
obj <- NormalizeData(obj) |>
  FindVariableFeatures(nfeatures = 200) |>
  ScaleData()

obj <- RunPCA(obj, verbose = FALSE) |>
  RunUMAP(dims = 1:10)

head(Embeddings(obj, "pca")[,1:3])
```

### 2.5 Graphs：KNN/SNN/WNN
```r
obj <- FindNeighbors(obj, dims = 1:10)
names(obj@graphs)  # 例如 "RNA_snn"
```

### 2.6 commands / tools / misc
```r
names(obj@commands)[1:5]           # 看看做過哪些主要步驟與參數
str(obj@tools, max.level = 1)
str(obj@misc,  max.level = 1)
```

### 2.7 小抄：我要拿 X，去哪裡？
| 需求           | 存放位置/類別     | 推薦取用方式                                  |
|----------------|--------------------|-----------------------------------------------|
| 原始計數       | `Assay@counts`     | `GetAssayData(obj, slot = "counts")`          |
| 正規化表達     | `Assay@data`       | `GetAssayData(obj, slot = "data")`            |
| 縮放後表達     | `Assay@scale.data` | `GetAssayData(obj, slot = "scale.data")`      |
| 高變異基因     | Assay/物件屬性     | `VariableFeatures(obj)`                        |
| 細胞中介資料   | `obj@meta.data`    | `obj$欄名` 或 `obj@meta.data`                 |
| 群集/身分      | Idents             | `Idents(obj)` / `RenameIdents()`              |
| PCA/UMAP 座標  | DimReduc           | `Embeddings(obj, "pca"/"umap")`               |
| KNN/SNN 圖     | `obj@graphs`       | `FindNeighbors()` 後 `names(obj@graphs)`      |
| 流程紀錄       | `obj@commands`     | `names(obj@commands)`                          |
| 細胞/基因名稱  | Cells/Features     | `Cells(obj)` / `Features(obj)`                |

---

## 3. 快速上手（pbmc_small）

目的：跑通最精簡流程；建立對 Seurat 物件的直覺。

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

## 4. 載入資料（10x/Matrix/RDS）

```r
# 10x 目錄（含 barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz）
# counts <- Read10X(data.dir = "data/10x_pbmc/")
# obj <- CreateSeuratObject(counts, project = "PBMC", min.cells = 3, min.features = 200)

# RDS（建議儲存/載入 Seurat 物件）
# saveRDS(obj, "data/seurat_obj.rds")
# obj <- readRDS("data/seurat_obj.rds")
```

---

## 5. 品質管制（QC）與過濾

```r
mt_pat <- if (any(grepl("^MT-", rownames(obj)))) "^MT-" else "^mt-"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pat)

VlnPlot(obj, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(obj, feature1 = "nCount_RNA",   feature2 = "nFeature_RNA") +
FeatureScatter(obj, feature1 = "nCount_RNA",   feature2 = "percent.mt")

# 以 MAD 為參考估計過濾閾值（請視資料調整）
mad_cut <- function(x, nmads=3){
  med <- median(x); m <- mad(x)
  c(max(0, med - nmads*m), med + nmads*m)
}
# lim_feat  <- mad_cut(obj$nFeature_RNA)
# lim_count <- mad_cut(obj$nCount_RNA)
# mt_upper  <- 20  # 人類常見 < 10–20%
# obj <- subset(obj, subset =
#   nFeature_RNA > lim_feat[1] & nFeature_RNA < lim_feat[2] &
#   nCount_RNA   > lim_count[1] & nCount_RNA   < lim_count[2] &
#   percent.mt   < mt_upper)
```

> 經驗閾值（依樣本差異調整）：`nFeature_RNA` 過低疑空滴、過高疑 doublet；人類 `percent.mt` 常見 < 10–20%。  
> 進階：`DropletUtils::emptyDrops` 去空滴；`SoupX` 校正環境 RNA。

---

## 6. 正規化與特徵選擇：SCTransform vs LogNormalize

**推薦：`SCTransform()`**（穩定、可回歸 `percent.mt`/細胞週期分數等）。
```r
# SCTransform（建議）
# obj <- SCTransform(obj, variable.features.n = 3000, vars.to.regress = "percent.mt", verbose = FALSE)
# DefaultAssay(obj) <- "SCT"
```

**傳統流程（適合教學/小資料）：**
```r
# LogNormalize
# obj <- NormalizeData(obj)
# obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
# obj <- ScaleData(obj, vars.to.regress = "percent.mt")
# DefaultAssay(obj) <- "RNA"
```

---

## 7. 降維、鄰居圖、分群與 UMAP

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

## 8. Marker 與初步註解（手動）

```r
# DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else "RNA"
# Idents(obj) <- "seurat_clusters"
# markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# head(markers)
# VlnPlot(obj, features = c("MS4A1","LYZ"), pt.size = 0)
# FeaturePlot(obj, features = c("MS4A1","CD3D","LYZ","GNLY"))
```

---

## 9. 自動註解（SingleR，可選）

```r
# library(SingleR); library(celldex); library(SingleCellExperiment)
# ref <- celldex::HumanPrimaryCellAtlasData()
# sce <- as.SingleCellExperiment(obj)
# pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
# obj$SingleR_label <- pred$labels
# DimPlot(obj, group.by = "SingleR_label", label = TRUE, repel = TRUE)
```

---

## 10. 雙重體偵測（scDblFinder）

```r
# library(scDblFinder); library(SingleCellExperiment)
# sce <- as.SingleCellExperiment(obj)
# sce <- scDblFinder(sce)
# table(sce$scDblFinder.class)
# obj$doublet <- sce$scDblFinder.class
# obj <- subset(obj, subset = doublet == "singlet")
```

---

## 11. 多樣本/批次整合（SCT / RPCA）

**SCT 整合（建議）**
```r
# features  <- SelectIntegrationFeatures(object.list = list_objs, nfeatures = 3000)
# list_objs <- PrepSCTIntegration(object.list = list_objs, anchor.features = features, verbose = FALSE)
# anchors   <- FindIntegrationAnchors(object.list = list_objs, normalization.method = "SCT",
#                                     anchor.features = features, verbose = FALSE)
# integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
# DefaultAssay(integrated) <- "integrated"
```

**RPCA 整合（適合相近組成/較大資料）**
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

## 12. 細胞週期與回歸

```r
# s.genes   <- Seurat::cc.genes.updated.2019$s.genes
# g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
# obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# obj <- SCTransform(obj, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = FALSE)
```

---

## 13. 差異表達（DE） 與 GSEA

**差異表達**
```r
# Idents(obj) <- "seurat_clusters"
# de_0_vs_all <- FindMarkers(obj, ident.1 = 0, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
# head(de_0_vs_all)
```

**GSEA（`clusterProfiler` + `fgsea` + `msigdbr`）**
```r
# library(clusterProfiler); library(msigdbr); library(fgsea)
# species <- "Homo sapiens"
# msig <- msigdbr(species = species, category = "H") |> dplyr::select(gs_name, gene_symbol)
# path_list <- split(msig$gene_symbol, msig$gs_name)
# mk <- FindMarkers(obj, ident.1 = 0, min.pct = 0.1, logfc.threshold = 0.1)
# rk <- mk$avg_log2FC; names(rk) <- rownames(mk); rk <- sort(rk, decreasing = TRUE)
# fg <- fgsea(pathways = path_list, stats = rk, minSize = 15, maxSize = 500, nperm = 1000)
```

> 常見錯誤：「`No gene can be mapped`」→ 檢查**物種**、**ID 類型**（`SYMBOL`/`ENTREZID`/`ENSEMBL`）與**大小寫**（人類多大寫，小鼠多小寫/首字大寫）。

---

## 14. 亞群重分析（Subsetting & Reclustering）

```r
# tcell_genes <- c("CD3D","CD3E","CD2")
# obj$TcellFlag <- Matrix::colSums(GetAssayData(obj, slot = "data")[tcell_genes, , drop = FALSE]) > 0
# tcell <- subset(obj, subset = TcellFlag)
# tcell <- SCTransform(tcell, variable.features.n = 2000, verbose = FALSE)
# tcell <- RunPCA(tcell); tcell <- FindNeighbors(tcell, dims = 1:20)
# tcell <- FindClusters(tcell, resolution = 0.6); tcell <- RunUMAP(tcell, dims = 1:20)
```

---

## 15. 資料導出與跨工具互通

```r
# saveRDS(obj, "export/seurat_obj.rds")

# library(SeuratDisk)
# SaveH5Seurat(obj, filename = "export/seurat_obj.h5seurat", overwrite = TRUE)
# Convert("export/seurat_obj.h5seurat", dest = "h5ad", overwrite = TRUE)  # for Scanpy
```

---

## 16. 可重現性與常見陷阱

```r
# if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
# renv::init(); renv::snapshot()
# sessionInfo()
```

**常見陷阱**
- 粒線體/核糖體前綴：Human `^MT-`/`^RP[LS]`；Mouse `^mt-`/`^Rp[ls]`。  
- 過濾過嚴/鬆：先看分佈再定閾值；MAD 只作參考。  
- 整合後仍有批次效應：調整 `nfeatures`/`dims`，或試 Harmony/fastMNN；先移除極端群。  
- 雙重體：提高預估率或調參重跑；檢查高 `nCount_RNA` 與多 lineage marker 的細胞。  
- GSEA 映射：物種/ID/大小寫一致；處理重複符號。  

---

## 附錄A：Minimal Pipeline 模板

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

## 附錄B：名詞對照（Glossary）

- **Feature/Gene**：基因（行）；**Cell/Barcode**：細胞（列）。  
- **Assay**：資料層（如 RNA、SCT、integrated）。  
- **Reduction**：降維結果（PCA/UMAP/TSNE）。  
- **Graph**：細胞間鄰居關係（KNN/SNN/WNN）。  
- **Anchors**：整合與標籤轉移的對應關係。  
- **Doublet**：兩個細胞共用同一條 barcode 的混合體。  

---

版權：本教學可自由使用與修改（CC BY 4.0）。
