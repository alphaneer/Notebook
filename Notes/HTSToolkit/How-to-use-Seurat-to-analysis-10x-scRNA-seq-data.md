# 如何使用Seurat分析10X Genomics的单细胞测序分析结果

10X Genomics数据经Cellranger流程定量后，下一步就推荐用R/Seurat进行探索性下游分析，这里用的数据来自于10X Genomics提供的外周血单核细胞(Peripheral Blood Mononuclear Cells, PBMC)的上游分析结果，共2,700个细胞被测序，测序平台为NextSeq 500,原始数据可通过<https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz>链接下载。

后续操作在R语言上进行，并假设工作路径和解压后的文件夹同级。

## 数据导入

10X Genomics分析流程CellRanger分析之后会有一个文件夹叫做"outs", 要导入的文件夹就在outs目录下。`Seurat`提供了`Read10X`专门导入这类数据。

```R
library(Seurat)
library(dplyr)
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, project = "10X_PBMC")
```

Seurat的分析需要构建`SeuratObject`对象，这一步可以对数据初步过滤，比如说只保留至少在3个细胞种表达的基因(min.cells)，只保留每个细胞至少要有200个基因被检测到。

## 标准预处理流程

### 质量控制和细胞筛选

在创建`SeuratObject`对象这一步已经做了初步的过滤，不过主要用途是节约内存，真正的过滤还需要根据探索性分析得到的QC图，才能选择更加合适的标准，从技术或生物学角度进行过滤。比如说对基因和分子量之间的关系进行可视化，剔除那些实验过程中造成的多个细胞被标记成单个barcode的情况。也可以基于线粒体基因比例来过滤细胞。

```R
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value =TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes,]) / Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")

```

基因和细胞的过滤需要根据具体试验进行，如下是一些文章用到的策略：

- Suijuan Zhong(2018 Nature): 仅保留基因表达量超过1000的细胞，仅保留标准化后表达量水平高于1，且至少在3个细胞中表达的基因。(从2394个细胞过滤到2333个细胞，留下17,854个基因)。之后进行初步聚类，剔除血红蛋白基因和小神经胶质特异的基因富集的两个聚类。次级聚类后，剔除了血红蛋白基因和免疫细胞。最后只有2309个细胞(不含免疫细胞)，且剔除了血红蛋白基因的结果用于t-SNE分析。
- Alexander(2018 Science): 剔除超过1% mt-RNA比对率的细胞。过滤表达基因数目低于250，且超过2500的细胞。初步聚类后，过滤在putative doublet clusters的细胞。层次聚类后，低于40个细胞的聚类剔除。
- Urszula(2018 eLife): 剔除没有表达，以及在对照组也表达的基因。 剔除低于25000 read counts以及低于1000基因表达的细胞。之后剔除在5个细胞中的reads小于10的基因。

### 标准化

### 检测单细胞间变化基因

### 数据缩放和剔除不需要的变异源