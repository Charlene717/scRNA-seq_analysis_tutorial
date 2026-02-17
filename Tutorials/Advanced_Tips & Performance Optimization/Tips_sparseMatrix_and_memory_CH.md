# sparseMatrix（稀疏矩陣）的優勢與使用情境（以單細胞/空間轉錄體為例）

`sparseMatrix`（R 中常見於 **Matrix** 套件，例如 `dgCMatrix`）最大的優勢是：**當資料中大部分元素是 0 時，可以大幅節省記憶體並提升運算效率**。  
這對 **scRNA-seq** 與 **spatial transcriptomics** 特別重要，因為 expression matrix 往往高度稀疏（常見 90–99% 為 0）。

<p align="center">
  <img src="https://raw.githubusercontent.com/Charlene717/scRNA-seq_analysis_tutorial/main/Figures/20260217_sparseMatrix_CH.png" width="900">
</p>
---

## 為什麼需要 sparse matrix？

以 scRNA-seq 為例，矩陣可能長這樣：

- 30,000 genes  
- 50,000 cells  
- 1.5 billion entries  
- 但真正有 expression 的可能只有 1–5%

如果使用一般 dense matrix（`matrix()`）：
- **所有 0 也會佔記憶體**

如果使用 sparse matrix（例如 `dgCMatrix`）：
- 只儲存：
  - 非 0 的數值
  - 這些數值的位置（row/column index）

---

## 核心優勢

### 1) 記憶體大幅節省
**Dense matrix** 需要儲存所有元素（包含 0），記憶體消耗非常可觀。  
**Sparse matrix** 只儲存非 0 值與索引，稀疏度越高，節省越明顯。

在單細胞與空間資料中，常見情況是：
- dense：數 GB～數十 GB
- sparse：數百 MB～約 1 GB（依資料規模而定）

---

### 2) 特定運算更快
當演算法能「跳過 0」時，sparse matrix 的計算量會顯著下降，例如：
- `rowSums()` / `colSums()`
- 矩陣乘法 `%*%`
- 部分 PCA / 線性代數操作（視實作而定）

**注意：**若某些步驟會把 sparse 強制轉成 dense，反而可能導致記憶體爆掉。

---

### 3) 非常適合高維資料
sparse matrix 特別適用於「高維 + 極度稀疏」的資料型態，例如：
- scRNA-seq
- ATAC-seq
- spatial transcriptomics（Visium HD / Xenium 等）
- graph adjacency matrix
- recommender system 類型的稀疏資料

---

## 在單細胞/空間分析中的實際好處

以 Seurat 為例，常見的 counts slot 本身就是 sparse：

- 可支援更大的細胞數（例如 100k+）
- 多樣本整合更可行
- 記憶體壓力顯著降低
- 許多常用步驟在稀疏表示下更有效率

---

## sparse matrix 的常見格式（dgCMatrix）

`dgCMatrix` 通常是「Compressed Sparse Column（CSC）」格式：
- 對「按 column 做運算」的效率很好（例如 `colSums()` 常很快）
- 也是多數單細胞框架常用的格式

---

## 什麼時候不適合用 sparse？

以下情境可能不適合或收益較小：
1. 矩陣本身不稀疏（例如 bulk RNA-seq 的 logCPM）
2. 大量 element-wise 操作且會頻繁改動矩陣結構
3. 某些套件或函式不支援 sparse，可能會暗中轉成 dense

---

## 常見陷阱（重要）

避免直接把 sparse 轉成 dense，例如：

- `as.matrix(sparse_obj)`

這可能瞬間把記憶體吃爆，尤其在單細胞與空間資料的實務規模下。

---

## 一句話總結
`sparseMatrix` 的核心優勢是：  
**在高度稀疏資料中，用更少的記憶體保存矩陣，並在許多運算上帶來更好的效率，是單細胞與空間轉錄體分析的基礎結構。**
