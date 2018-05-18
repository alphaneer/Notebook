# RNA-Seq工作流：基因水平探索性分析和差异表达

----
[TOC]

## 实验数据

数据存放在`airway`包中，主要是气道平滑肌细胞在DEX处理后。糖皮质激素被哮喘患者用于减少气道炎症。
四组处理和四组未处理。
文献为： [PubMed entry 24926665](http://www.ncbi.nlm.nih.gov/pubmed/24926665)
GEO entry ：[GSE52778](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778).

## count matrix准备

count-based数理统计包要求得到数据格式为，第i列j行为多少read。
数据**不能**根据测序深度和文库预先进行**标准化**

### 新的方法

**注意**：目前有转录本丰度定量方法，如[Salmon](https://combine-lab.github.io/salmon/)(Patro et al. 2016), [Sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/) (Patro, Mount, and Kingsford 2014), [kallisto](https://pachterlab.github.io/kallisto/) (Bray et al. 2016), or [RSEM](http://deweylab.github.io/RSEM/) (B. Li and Dewey 2011),可以直接估计丰度，而不必进行序列比对，之后使用`tximport`包用于组装序列然后生成差异化表达分析的矩阵。

优点：

- 更加高效
- 避免了基于比对会导致基因比对到多个同源基因的情况

### 将读段比对到参考基因组

根据benchmarking paper找合适的比对软件，此处是用STAR

```bash
for f in `cat files`; do STAR --genomeDir ../STAR/ENSEMBL.homo_sapiens.release-75 \
--readFilesIn fastq/$f\_1.fastq fastq/$f\_2.fastq \
--runThreadN 12 --outFileNamePrefix aligned/$f.; done
```

SAMtools用于转换格式，

```bash
for f in `cat files`; do samtools view -bS aligned/$f.Aligned.out.sam \
-o aligned/$f.bam; done
```

### 练习所用数据

实际操作为了方便只选取了部分数据，

```r
library("airway")
dir  <-  system.file("extdata",package  =  "airway",mustWork  =  TRUE)
list.files(dir)
```

每个项目需要为每一个样品提供详细的信息，如

```r
csvfile <- file.path(dir,"sample_table.csv")
(sampleTable <- read.csv(csvfile,row.names = 1))
```

**小技巧**：在赋值操作外加括号，可以在赋值操作结束后打印结果到屏幕。

### 根据比对文件统计read/fragment数量：

#### 可用工具：

| function | package | framework | output | _DESeq2_ input function |
| :-- | :-- | :-- | :-- | :-- |
| _summarizeOverlaps_ | _[GenomicAlignments](http://bioconductor.org/packages/GenomicAlignments)_ | R/Bioconductor | _SummarizedExperiment_ | _DESeqDataSet_ |
| _featureCounts_ | _[Rsubread](http://bioconductor.org/packages/Rsubread)_ | R/Bioconductor | matrix | _DESeqDataSetFromMatrix_ |
| _tximport_ | _[tximport](http://bioconductor.org/packages/tximport)_ | R/Bioconductor | list of matrices | _DESeqDataSetFromTximport_ |
| _htseq-count_ | [HTSeq](http://www-huber.embl.de/users/anders/HTSeq) | Python | files | _DESeqDataSetFromHTSeq_ |

```r
# 提供文件位置
filenames <- file.path(dir,paste0(sampleTable$Run,"_subset.bam"))
file.exists(filenames)
# 导入
library("Rsamtools")
bamfiles <- BamFileList(filenames,yieldSize = 2000000)
seqinfo(bamfiles[1])
```

**注意**：请保证注释gtf文件里的染色体命名和比对得到文件染色体命名一致。

#### 提供基因模型

```r
library(GenomicFeatures)
gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")
(txdb <- makeTxDbFromGFF(gtffile,format="gtf",circ_seqs = character()))
# 这一步可用产生所有外显子的范围RangeList
(ebg <- exonsBy(txdb, by="gene"))
```

#### Read counting

使用GenomicAlignments的summarizeOverlaps对bamfiles中read进行count。

```r
library("GenomicAlignments")
library("BiocParallel")
register(SerialParam())
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
```

**注意**：是否单端还是双端，文库有没有正负链特异性

#### SummarizedExperiment对象

![SummarizedExperiment对象](http://www.bioconductor.org/help/workflows/rnaseqGene/rnaseqGene_files/figure-markdown_strict/sumexp-1.png)

assay存储counts矩阵,rowRanges是基因组范围(genomic ranges)信息,colData则是样品信息, 高亮的部分为第一行

以我们的se为例：

```r
class(se)
dim(se)
assayNames(se)
head(assay(se),3)
rowRanges(se)
str(metadata(rowRanges(se)))
colData(se)

(colData(se) <- DataFrame(sampleTable))
```

### 分支点

获得count统计数据后的分析有许多包可以用来后续分析：

- _[edgeR](http://bioconductor.org/packages/edgeR)_ (M. D. Robinson, McCarthy, and Smyth 2009),
- _[limma](http://bioconductor.org/packages/limma)_ with the voom method (Law et al. 2014),
- _[DSS](http://bioconductor.org/packages/DSS)_ (H. Wu, Wang, and Wu 2013),
- _[EBSeq](http://bioconductor.org/packages/EBSeq)_ (Leng et al. 2013)
- _[baySeq](http://bioconductor.org/packages/baySeq)_ (Hardcastle and Kelly 2010)

Schurch比较不同RNA-Seq统计方法的表现，可以通过[这篇文章](https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/27022035/)了解如何进行实验设计和工具选择

## DESeqDataSet：样品信息和设计矩阵(design matrix)

DESeqDataSet基于_SummarizedExperiment_，有几点不同

- the `assay` slot is instead accessed using the _counts_ accessor function, and the _DESeqDataSet_ class enforces that the values in this matrix are non-negative integers.
- _DESeqDataSet_ has an associated _design formula_.

design matrix也就是实验设计，他为DESeq2函数在分析时如何处理样本提供信息。
差异表达分析所需要的design formula存储在colData，colData需要提供表型数据

### 从SummarizedExperiment开始

此处采用使用数据产生的结果

```r
data("airway")
se <- airway
# 改变factor的排序
se$dex <- relevel(se$dex,"untrt")
# 查看design matrix信息
colData(se)
dds <- DESeqDataSet(se, design = ~ cell + dex)
```

### 从count matrix开始

```r
countdata <- assay(se)
head(countdata,3)
coldata <- colData(se)
ddMat <- DESeqDataSetFromMatrix(countData=countdata,colData=coldata,design= ~cell+dex)
```

## 探索性分析和可视化

### 数据集预处理

删除那些数据量过少

```r
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1,]
nrow(dds)
```

### rlog（regularized-logarithm transformation）变换

许多常见的多维数据探索性分析的统计分析方法，例如聚类和主成分分析，能在某一类数据中有更好的效果，这类数据虽然平均值不同，但是方差相同。当方差期望值在不同均值下几乎相同时，数据可以称为**等方差性**。然而，对于RNA-Seq原始统计，反差伴随均值增加。如果直接在尺寸因子标准化的read统计上进行PCA分析，结果通常会依赖于少数几个高表达基因，因为他们在不同样本间具有最大的绝对差。一个比较简单常用的方法是在标准化的count值得对数上加1个1个伪count。尽管，取决于选择什么样的伪count，但是现在那些count非常低的基因对于结果作图也具有相同的噪点，因为低countd的对数实际上扩大了他们的方差。我们可以通过一些拟真数据（Poisson counts，lambda从0.1到100）可视化这种特性。

以平均数为x轴做每行（不同基因）的标准离差分布图

```r
lambda <- 10^(seq(from=-1,to=2,length=1000))
cts <- matrix(rpois(1000*100,lambda),ncol=100)
library("vsn")
library("hexbin")
meanSdPlot(cts,ranks=FALSE)
# 变换后
log.cts.one <- log2(cts+1)
meanSdPlot(cts,ranks=FALSE)
```

DESeq2为count数据提供了两类变换方法，使得不同均值的方差趋于稳定：**regularized-logarithm transformation or rlog**和**variance stabilizing transformation(VST)**,含有色散平均趋势负二项数据。

对于高count的基因，rlog和VST能提供相似的结果。对于低count的基因，值会根据不同样本基因的平均值发生收缩。之后VST和rlog变换后数据差不多是等方差了，就可以直接用于计算不同样本间的距离，PCA作图或者用作下游分析的输入。

**方法选择：**rlog,数据集小于30，VST大数据集。

**两类变换不能用于差异检验**，差异检验需要使用原始数据

```r
# rlog变换
rld <- rlog(dds,blind=FALSE)
head(assay(rld),3)
# vst变换
vsd <- cst(dds,blind=FALSE)
head(assay(vsd),3)
```

blind=FALSE,意为cell line和treatment之间的差异不会作用于实验期望的方差-平均值趋势。

```r
par(mfrow=c(1,3))
dds <- estimateSizeFactors(dds)
lims <-c(-2,20)
plot(log2(counts(dds,normalized=TRUE)[,1:2]+1),pch=16,cex=0.3,main="log2(x+1)",xlim=lims,ylim=lims)
plot(assay(rld)[,1:2],,pch=16,cex=0.3,main="rlog",xlim=lims,ylim=lims)
plot(assay(vsd)[,1:2],,pch=16,cex=0.3,main="vst",xlim=lims,ylim=lims)
```

![双样本变换后count的散点图](http://www.bioconductor.org/help/workflows/rnaseqGene/rnaseqGene_files/figure-markdown_strict/rldplot-1.png)

### 样本距离

RNA-Seq分析**第一步**通常是评估样本间的总体相似度。

- 那些样本最接近
- 那些样本差异较大
- 这与实验设计预期符合么
- 使用R内置的 _dist_ 计算 Euclidean distance

```r
sampleDists <- dist( t( assay(rld) ) )
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

- _Poisson Distance_, 泊松距离同样考虑到counts的内在固有方差结构

```r
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(samplePoisDistMatrix) <- NULL

pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)
```

#### PCA plot

还有一种可视化样本-样本距离的方法就是主成分分析。在这个排序方法中，数据点（样本）投影到2维平面上，因此他们在2个方向上分散情况解释了大多数的差异，

```r
## 使用DESEeq2的plotPCA
plotPCA(rld, intgroup = c("dex", "cell"))

## 使用ggplot2
library("ggplot2")
pcaData <- plotPCA(rld, intgroup = c( "dex", "cell"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=dex, shape=cell)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
```

#### MDS plot

当我们没有原始数据数组时，但是有距离数组，可以使用MDS（multidimensional scaling)

```r
# rlog变换counts计算距离
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds,aes(X1,X2,color=dex,shape=cell)) + geom_point(size=3) + coord_fixed()
```

泊松距离用MDS作图

```r
mdsPoisData <- data.frame(cmdscale(samplePoisDistMatrix))
mdsPois <- cbind(mdsPoisData, as.data.frame(colData(dds)))
ggplot(mdsPois, aes(X1,X2,color=dex,shape=cell)) + geom_point(size=3) +
  coord_fixed()
```

## 差异表达分析

### 差异分析管道

由于我们已经声明了实验设计在创建_DESeqDataSet_时，我们可以在原始counts运行差异表达pipeline，只需要调用函数DESeq

```r
dds <- DESeq(dds)
```

具体步骤见?DESeq：estimation of size factors,the estimation of dispersion values for each gene, fitting a generalized linear modle

_DESeqDataSet_返回内容含有所有匹配的参数。

### 构建结果表格

不加其他参数调用_results_会提取预测log2倍变化和p值

```r
res <- results(dds)
```

res为DataFrame对象，带有各列含义信息的元数据

```r
mcols(res, use.names=TRUE)

## DataFrame with 6 rows and 2 columns
##                        type                               description
##                 <character>                               <character>
## baseMean       intermediate mean of normalized counts for all samples
## log2FoldChange      results  log2 fold change (MAP): dex trt vs untrt
## lfcSE               results          standard error: dex trt vs untrt
## stat                results          Wald statistic: dex trt vs untrt
## pvalue              results       Wald test p-value: dex trt vs untrt
## padj                results                      BH adjusted p-values
```

baseMean：标准化count值的平均值，除以size factors，所有DESeqDataSet的样本。随后的四列为特定的比较。如trt水平与untrt水平在dex上的比较。我们在后面会支出如何得到其他比较。

log2FoldChange: 效应大小估计，effect size estimate.基因表达多大程度是不同dex处理引起的。

对结果描述性分析后，会提供额外的信息。`summary(res)`

许多因为dex处理而差异表达的基因的FDR（假阳性率）在10%。这很正常，因为气管平滑肌细胞已知能对糖皮质激素类固醇产生反应。有两种方法可以限制显著性的基因

- 降低FDR阈值
- 提高log2 fold change阈值，ifcThreshold

如果我们要降低了FDR阈值，我们需要将该值告诉results()

```r
result.05 <- results(dds,alpha=0.05)
table(res.0.05$padj) <.05)
```

如果要提高log2 fold change threshold，这样能发现那些基因是的确因为不同处理而发生变化。

```r
resLFC1 <- results(dds,lfcThreshold=1)
table(resLFC1$padj <0.1)
```

citation("pkgName")可以知道如何对包进行引用。

#### 其他比较

一般而言，任意两个水平的变量比较结果可以用_contrast_参数取得。使用者需要声明3个值：变量名、分子水平的名字，分母水平的名字

```r
results(dds, contrast=c("cell","N061011","N61311"))
```

#### 多重测试

在高通量生物学，我们会比较谨慎的使用p值否定null假设，但是用于矫正多重测试。当我们把p值调到0.05，就只有5648个基因

```r
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))
```

假设所有基因的零假设是真实的，即，没有基因受到地塞米松的治疗的影响。那么通过定义p值，我们预期5%的基因p值低于0.05是差异表达的，那么久有1470/5648 =25%的假阳性。
DESeq2使用Benjamini-Hochberg (BH) 修正在基础Rp.adjust函数中使用。
FDR在许多高通量试验中是一个非常有用的统计值，因为我们经常需要报告或注意一组感兴趣基因，我们会在这组数据上加一个FDR上限
因此，如果我们考虑到10%的假阳性是可接受的，我们就可以认为所有修正p值小于10%=0.1的基因是显著的，那么有多少基因呢？

```r
sum(res$padj<0.1,na.rm=TRUE)
```

提取这些基因并根据log2 fold change进行排序 ，找到明显下调显著性的基因

```r
resSig <- subset(res,padj<0.1)
head(resSig[order(resSig$log2FoldChange),])
```

明显上调的基因

```r
head(resSig[ order(resSig$log2FoldChange,decrease=TRUE),])
```

## 结果作图

### 单个基因不同处理组的标准化counts

快速对特定基因counts进行可视化的函数是_plotCounts_，它接受DESeqDataSet，基因名，以及用于作图的组作为参数

```r
topGene <- rownames(res)[which.min(respadj)]
plotCounts(dds,gene=topGene,intgroup=c("dex))
```

### 用颜色声明细胞系

使用ggplot2的ggplot函数自定义画图

```r
geneCounts <- plotCounts(dds,gene=topGene,intgroup=c("dex","cell"),returnData=TRUE)
ggplot(geneCounts,aes(x=dex,y=count,color=cell)) + scale_y_log10() + geom_point(position=position_jitter(width=.1,height=0),size=3)
```

### 使用更加结构化的重排

```r
ggplot(geneCounts, ae(x=dex,y=count,fill=dex)) + scale_y_log10() + geom_dotplot(binaxis="y",stackdir="center")
```

### Normalized counts with lines connecting cell lines

```r
ggplot(geneCounts, aes(x=dex, y=count, color=cell, group=cell)) +
  scale_y_log10() + geom_point(size=3) + geom_line()
```

### An MA-plot of changes induced by treatment

MA-plot为一个两组比较试验提供了好用的总览.
This plot demostrates that only genes /with a large average normalized count /contain sufficient information /to yield a significant call.

```r
plotMA(res,ylim=c(-5,5))
```

### p值柱形图

平均标准化count值大于1的基因的p值柱状图

```r
hist(res$pvalue[res$baseMean > 1 ], breaks=0:20/20, col="grey50",border="white")
```

### 基因聚类

之前的树状图提供了样本之前的层级关系。这种层级关系也能应用于基因水平。Since the clustering is only relevant for genes that actually carry a signal, one usually would only cluster a subset of the most highly variable genes.

```r
library(genefilter)  # 载入genefilter用于筛选
topVarGenes <- head(order(rowVars(assay(rld)),decreasing = TRUE),20)
# rowVars{genefilter}:数值型数组行方差
# 也就是找出前20个不同样本（经过标准化）间方差，降序，返回行号
mat <- assay(rld)[ topVarGenes, ]
# 将具体的样本的标准化count传递给mat
mat <- mat - rowMeans(mat)
# 减去平均值
df <- as.data.frame(colData(rld)[,c("cell","dex")])
# 提供注释信息
pheatmap(mat,annotation_col  = df)
```

### independent filtering 独立筛选

MA plot强调了RNA-Seq数据一个重要特性。对于弱表达基因，我们没有机会观察差异表达，因为low read counts由于有太高的泊松噪点以至于任何生物学效应会以row rate淹没在不确定中。
we can also show this/ by examining the ratio of samll p values(say, less than 0.05)/ for genes binned by mean normalized count.
下列代码，首先使用_quantile_函数创建分组，使用_cut_根据base mean对基因进行分组，使用中位点重命名levels，为每一个分组计算p值低于0.05的比率，然后做图

```r
resLFC1 <- results(dds,lfcThreshold=1)
qs <- c(0,quantile(resLFC1$baseMean[resLFC1$baseMean >0],0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~",round(signif(.5*qs[-1]+.5*qs[-length(qs)],2))
ratios <- tapply(resLFC1$pvalue,bins,function(p) mean(p < .05,na.rm=TRUE))
barplot(ratios,xlab="mean normalized count", ylab="ratio of small p values")
```

在FDR过程中将输入移除low count genes，我们可以在我们所保留的数据中找到更多的显著性基因。

DESeq2基因能够自动进行independent filtering，最大化那些修正p值小于重要值的数量。可以在_result_函数进行控制。

## 注释和输出结果

### 注释

结果表目前只含有Ensembl gene ID信息，但是其他基因名对于合作者可能更具有信息性。Bioconductor's 注释包帮助我们让多个ID命名之间能够相互映射

```r
library("AnnotationDbi")
library("org.Hs.eg.db")
```

- org: organsim annotation
- Hs: Homo sapiens
- db: 以AnnotationDbi数据库方式组织
- eg: 使用Entrez Gene ID 作为主键

```r
# 可用的关键字类型
columns(org.Hs.eg.db)
```

可用_mapIds_添加独立列到结果表中，

```r
res$symbol <- mapIds(org.Hs.eg.db,keys = row.names(res),column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res$entrez <- mapIds(org.Hs.eg.db,keys = row.names(res),column="ENTREZID",keytype = "ENSEMBL",multiVals = "first")
```

这样结果就有了合适的外部基因ID

```r
resOrdered <- res[order(res$padj),]
head(resOrdered)
```

### 结果输出

将结果保存为csv文件

```r
resOrderedDF <- as.data.frame(resOrdered[1:100],)
write.csv(resOrderedDF,file="result.csv")
```

更加高级的工具为Reporting Tools。

```r
library("ReportingTools")
htmlRep <- HTMLReport(shortName="report", title="My report",
                      reportDirectory="./report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)
```

### 在基因组上画出fold change

如果使用_summarizedOverlaps_函数对reads技术，那么DESeqDataSet就已经为每个基因表明了range

```r
resGR <- results(dds,lfcThreshold=1,format="GRanges")
resGR
resGR$symbol <- mapIds(org.Hs.eg.db,names(resGR),"SYMBOL","ENSEMBL")
library("Gviz")

window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)

sig <- factor(ifelse(resGRsub$padj < .1 & !is.na(resGRsub$padj),"sig","notsig"))

```

```r
options(ucscChromosomeNames=FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name="gene ranges", feature=sig)
d <- DataTrack(resGRsub, data="log2FoldChange", baseline=0,
               type="h", name="log2 fold change", strand="+")
plotTracks(list(g,d,a), groupAnnotation="group", notsig="grey", sig="hotpink")
```