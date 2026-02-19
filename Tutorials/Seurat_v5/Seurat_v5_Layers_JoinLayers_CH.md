# JoinLayers 函數全面說明（Seurat v5）

## 1 背景：Seurat v5 的層（layers）概念

Seurat v5 引入了 **Assay5** 類別，允許在同一個 assay 內存儲多個 **層**。這些層通常用來表示不同樣本或處理條件的計數、正規化和縮放矩陣，例如在多樣本合併後，`RNA` assay 可能同時包含 `counts.sample1`、`counts.sample2`、`data.sample1` 等。這種設計可讓使用者在單個物件內保留樣本特異的資訊而不立即縫合矩陣。然而，舊版 Seurat（v4 及以前）採用的是單一矩陣結構，只有 `counts` 和 `data` 插槽。在 v5 中直接呼叫 v4 時代的函式（如 `GetAssayData`、`FindMarkers`）會因未合併層而報錯，這也是許多使用者誤以為 `JoinLayers()` 用於「適配舊版」的原因【245596585732741†L220-L233】。

## 2 函式概覽

### 2.1 主要作用

`JoinLayers()` 用於**將 Assay5 或 Seurat 物件中的多個層重新縫合成單一矩陣**。它會按照細胞名稱和基因名稱對齊樣本矩陣，再將它們串接起來，最終只保留 `counts`、`data` 及 `scale.data` 層（或使用者指定的層名稱）。這個動作不會執行任何批次校正或整合，也不會影響已經儲存在 reduction 槽的整合結果【352989248222699†L244-L252】。

### 2.2 語法與參數

對於 `Assay5`，`JoinLayers()` 的常見使用方式為：

```R
JoinLayers(object, layers = NULL, new = NULL, default = TRUE)
```

- **`object`**：要處理的 Assay5 或 Seurat 物件。若傳入 Seurat 物件，可透過 `assay` 參數指定要處理的 assay。
- **`layers`**：要合併的層名向量。預設值為 `counts`、`data`、`scale.data`【430440023730044†L680-L688】。
- **`new`**：合併後新層的名稱；若為 `NULL`，則使用舊名稱【430440023730044†L680-L688】。
- **`assay`**：僅在處理 Seurat 物件時使用，用以指定要合併的 assay【471935816826373†L18-L23】。
- **返回值**：合併後的 Assay5 或 Seurat 物件。

呼叫後，函式會為每個指定層調用 `JoinSingleLayers()` 進行縫合，並移除原層。如果 `default` 為 `TRUE`，新的層會成為 assay 的默認層【430440023730044†L1243-L1275】。

## 3 使用情境

### 3.1 合併或拆分樣本後

* **使用 `merge()` 後**：在 v5 中，`merge()` 會保留各樣本的 layer 而不是立即拼接矩陣【692994102636461†L254-L263】。若未計畫進行批次校正，可直接對合併後的物件執行 `JoinLayers()` 使其回到單層結構，例如：

  ```R
  merged_obj <- merge(x = sample1, y = list(sample2, sample3))
  merged_obj[['RNA']] <- JoinLayers(merged_obj[['RNA']])
  ```

* **使用 `split()` 或 `SplitObject()` 後**：將一個對象按照某個 metadata 拆分成多個 layer 以評估整合品質。完成後可用 `JoinLayers()` 將其重新合併【692994102636461†L234-L243】。

### 3.2 整合與批次校正流程中

Seurat v5 的整合（例如 `IntegrateLayers()` 或 `IntegrateData()`）在低維度空間計算，不會直接修改基因計數矩陣。整合結束後，為了進行差異表達分析，需要重新縫合 count 和 data 層。官方整合教程指出：「整合完成後，你可以重新合併 layers —— 這會將單獨的資料集合併，並重建原始的 `counts` 和 `data` 層。**在進行差異表達分析之前必須這樣做**」【438342738653614†L320-L324】。單細胞研討課程也提醒，整合後如要檢視群集間的差異，必須執行 `JoinLayers()`【245596585732741†L220-L233】。

### 3.3 與第三方函式或舊版工作流程相容

許多 v4 時期開發的函式（例如 `GetAssayData`、`FindMarkers`、`FindConservedMarkers`、SingleR、CellChat 等）預期 assay 僅有單一矩陣。若在 v5 中未合併層，這些函式會報錯「data layers are not joined」。SeuratExtend FAQ 建議在 `merge()` 合併樣本後立即執行 `JoinLayers()`，才能正常使用這些函式【801390038896824†L124-L146】。因此，`JoinLayers()` 是將 v5 對象格式調整為 v4 期望的結構的關鍵步驟，卻並不真正改變版本。

### 3.4 多模態資料

在 CITE-seq 等多模態資料中，一個 Seurat 對象可能含有 RNA、ADT 等多個 assay，每個 assay 內又有多個層。使用 `JoinLayers()` 時需要為每個 assay 分別處理，例如先對 `RNA` 執行，再對 `ADT` 执行【844843906040652†L38-L44】。

### 3.5 記憶體與性能

當對象包含大量樣本或多層時，層資料會佔用額外的記憶體。中國的解析文章指出，合併層可以減少 3–5% 的內存占用並簡化資料結構【841436903913337†screenshot】。此外，GitHub 用戶建議在合併大量對象前先對每個對象執行 `JoinLayers()`，這樣再 `merge()` 所有對象時效率更高【892034125168339†L389-L392】。

### 3.6 將 Seurat 對象轉換為其他格式

在將 Seurat v5 對象導出為 SingleCellExperiment (`sce`) 或 `.h5ad` 文件時，其他框架並不支持多層結構，需先透過 `JoinLayers()` 合併層，保證 `counts` 和 `data` 只有一個矩陣【245596585732741†L220-L233】。

## 4 與整合流程的關係

重要的是理解 `JoinLayers()` 並不會替代整合操作：

* **非整合函式**：開發者在 GitHub 說明，`JoinLayers` 只負責縫合多個樣本的計數矩陣，並不會執行任何批次校正或整合，也不會改變降維槽中的整合結果【352989248222699†L244-L252】。

* **整合後必需**：整合或 batch correction 後對象仍保留每個樣本的 layer。若直接對未合併的層執行差異表達分析，會引發錯誤；故必須先用 `JoinLayers()` 將層合併【438342738653614†L320-L324】【758681675104518†L343-L347】。

## 5 與舊版 Seurat 的相容性

`JoinLayers()` 常被誤認為是用於舊版相容，但它的設計重點在於處理 v5 的層結構：

* **不是降級工具**：`JoinLayers()` 只是將多個層串接成單一矩陣，以適配需要單層資料的函式和包。若想真正將 v5 物件「降級」為 v4/v3 格式，必須先執行 `JoinLayers()`，再透過 `as(object, Class = 'Assay')` 或 `CreateAssayObject()` 將 `Assay5` 轉成舊版 `Assay` 類別【626575938515563†L397-L417】。

* **保存 metadata**：當使用 `JoinLayers()` 合併層時，應將返回的 `Assay5` 重新賦值回原對象的 assay 插槽 (`seurat[['RNA']] <- JoinLayers(seurat[['RNA']])`)；若直接寫成 `seurat <- JoinLayers(seurat[['RNA']])` 會丟失 meta.data【367101772609323†L341-L347】。

* **SCTAssay 不適用**：如果默認 assay 是 `SCT`（SCTransform），無法直接對其使用 `JoinLayers()`。應先將 default assay 切換回 `RNA`，對 `RNA` assay 執行 `JoinLayers()`，再繼續分析【305772341404596†L24-L31】。

## 6 注意事項與進階建議

1. **核對層結構**：執行 `Layers(object)` 可查看當前 assay 中有哪些層，確認合併前是否存在多個樣本層。

2. **合併順序**：建議在所有樣本質控、初步合併後立即執行 `JoinLayers()`；若需要進行整合，可等整合完成後再合併層。下游差異分析必須在合併層後進行【758681675104518†L343-L347】。

3. **多模態合併**：對於含多個 assay 的對象，應分別對每個 assay 執行 `JoinLayers()`，以免忽略其他模態【844843906040652†L38-L44】。

4. **可逆性**：合併層後仍可透過 `split()` 或 `SplitObject()` 再次拆分層，以進行樣本特異分析【692994102636461†L234-L243】。

5. **記憶體管理**：若資料集龐大，建議先對每個對象執行 `JoinLayers()`，再進行大型 `merge()` 以減少執行時間【892034125168339†L389-L392】。此外，合併層可降低計算過程的記憶體消耗【841436903913337†screenshot】。

## 7 範例工作流程

以下示範從載入資料、合併樣本、整合到差異分析的完整流程：

```R
# 讀取兩個樣本
sample1 <- CreateSeuratObject(counts = Read10X("/path/to/sample1"), project = "S1")
sample2 <- CreateSeuratObject(counts = Read10X("/path/to/sample2"), project = "S2")

# 合併樣本（v5 保留多個 layer）
merged <- merge(sample1, y = list(sample2))

# （可選）執行批次校正：
# merged <- IntegrateLayers(merged, origidents = "project")

# 合併層（必須在差異表達分析前完成）
merged[['RNA']] <- JoinLayers(merged[['RNA']])

# 正規化、變異基因尋找與降維
merged <- NormalizeData(merged) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

# 聚類
merged <- FindNeighbors(merged, dims = 1:10) %>% FindClusters(resolution = 0.5)

# 差異表達分析
markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25)

# 進階：再次拆分層進行樣本特異分析
# merged_list <- split(merged, merged$project)
```

在整合流程中，應先質控並合併樣本（保留多層），然後進行整合 (例如使用 `IntegrateLayers()`)，評估整合效果後，再使用 `JoinLayers()` 合併層，最後進行差異分析和其他下游操作【438342738653614†L320-L324】。

## 8 總結

* **核心功能**：`JoinLayers()` 將 Seurat v5 中拆分或合併後的多個樣本層縫合為單一矩陣，使對象回到類似 v4 的結構。【758681675104518†L343-L347】強調它是批次校正和整合後進行差異分析的前置步驟。

* **用途多樣**：除了差異分析，它還能提升與舊版或第三方函式的相容性【801390038896824†L124-L146】，降低記憶體消耗【841436903913337†screenshot】，並在將資料導出到其他框架前確保格式正確【245596585732741†L220-L233】。

* **非整合工具**：`JoinLayers` 不會改變整合結果或執行批次校正；整合須由 `IntegrateLayers()` 等函式完成【352989248222699†L244-L252】。

* **注意保存 meta data**：務必將返回的 `Assay5` 重新放回原來的 assay 插槽，以免刪除 meta.data【367101772609323†L341-L347】。

透過理解 `JoinLayers()` 的角色與使用時機，研究者可以靈活管理 v5 的多層資料結構，在需要時合併或拆分層，既享有彈性又兼顧下游分析的要求。