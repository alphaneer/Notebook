# RNA-seq analysis is easy as 1-2-3 with limma, Glimma, and edgeR

---
[TOC]

# 摘要

Bioconductor项目的一个主要优势是简单和有效的对RNA-seq数据进行分析。从基因水平的counts summarised开始，典型的分析包括：预处理，探索性数据分析，差异表达检验和利用之后试验和验证研究得到的结果进行通路分析。在这篇工作流程文章中，我们分析了小鼠乳腺的RNA-Seq数据，使用广受欢迎的**edgeR**包对数据导入、组织、过滤和标准化，然后使用**limma**的_voom_工具的线性建模、经验贝叶斯修正评估差异表达并进行gene set检验.这个工作流程随后被**Glimma**所优化，它能交互探索结果，这样使用者能够验证单个样本和基因。这三个包所提供的完整分析流程强调了方便性，这样研究者就能使用Biconductor从生物学角度观察RNA-Seq试验得到的原始counts

# 简介

RNA-sequencing（RNA-Seq）目前已成为基因表达谱研究的主要技术，它能通过全基因组检测找到两个或多个条件下差异表达基因。Bioconductor项目的**edgeR** (Robinson, McCarthy, and Smyth 2010) 和 **limma** packages (Ritchie et al. 2015)为处理这些RNA-Seq数据问题提供了一套成熟的统计学方法。

在这篇文章中，我们描述了用于分析RNA-Seq数据的edgeR-limma工作路程，它使用基因水平的counts作为输入，随后通过预处理和探索性数据分析获得一系列差异表达（DE）基因和基因指纹。这个分析流程随后通过**Glimma**的交互图形界面得到提高，它允许对同时在样本和基因水平中对数据进行更加细致的套索，而不是仅仅能用R的统计作用。

该工作流程所分析的试验来自于Sheridan _et al._ (2015) (Sheridan et al. 2015)，包括三个细胞群体(basal, luminal progenitor (LP) and mature luminal (ML))，根据处女小鼠的乳腺排序，每一个都有三个重复。RNA样本在Illumina HiSeq 2000分为三个批次获得100bp单端(PE)reads(读段)。本文的分析要点认为RNA-Seq试验得到的reads已经比对到恰当的比对基因组并且概要到基因特异性区域关联counts。这个实例中，reads已经比对到小鼠参考基因组(mm10)，使用**R**基本流程的**Rsuberead**包（align函数后用featureCounts,for gene-level summarisation based on the in-built mm10 RefSeq-based annotation)

样本的Count数据可以从Gene Expression Ommnibus(GEO)下载，序列号为GSE63310.更多实验设计和样本准备信息见GEO相应的序列号下。

# 配置

```r
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
```

# 数据包装

## 读取count-data

下载GSE63310_RAW.tar并解压

```r
url <- "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
utils::untar("GSE63310_RAW.tar", exdir = ".")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
  "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
  "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)
```

每一个文本文件都包含给定样本的原始基因水平counts。我们的分析只需要试验中basal,LPhe ML样本。

```r
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", 
   "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt", 
   "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", 
   "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", 
   "GSM1545545_JMS9-P8c.txt")
read.delim(files[1], nrow=5)
```

尽管这9个文本文件可以分别读入R中然后再合成一个counts矩阵，**edgeR**有一种更方便的方法，readDEG.结果会生成一个DGEList对象，包括27,179行和对应的Entrez gene ID,9列试验的样本个体。

```r
x <- readDGE(files, columns=c(1,3))
dim(x)
class(x)
```

如果count数据在一个文件中，则用`DGElist`函数

## 组织样本信息

> 光有count matrix还不够，你还需要补充实验设计的一些信息，比如说试验变量

为了下游分析，试验设计有关的样本水平信息需要关联到count matrix的列中。这应该包括对表达水平有影响的试验变量（生物学和技术上）。样本包括细胞类型(basal, LP和 ML),基因型（wild-type,knock-out)，表现型(疾病状态，sex,age)样本处理(drgu,control)和批次信息（如果样本在不同时间点上手气和分析，那么实验日期需要注明），以上仅是一部分。

```r
DGEList(counts = matrix(0, 0, 0), lib.size = colSums(counts),
        norm.factors = rep(1,ncol(counts)), samples = NULL,
        group = NULL, genes = NULL, remove.zeros = FALSE)
```

系统默认生成的DESList对象中的sample，是一个data.frame,可以添加试验设计的有关数据

我们得到DESList-object包含一个`samples`数据框，保存cell type和batch（测序泳道）信息，每一个都有三个不同水平。`x$samples`的文库大小会自动计算，标准化因子会设置为1.为了简明起见，我们修改了一下DGEList-object x的列名。

```r
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames

colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP",
                     "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
x$samples
```

> 上述操作就是在x的samples修改了样品的命名，添加了分组信息和泳道信息

## 组织基因注释

DGEList-object第二个data.frame `genes`用于储存gene-level信息，与counts matrix的行对应。这类信息来源有
1.organism specific packages such as **Mus.musculus** (Bioconductor Core Team 2016b) for mouse (or **Homo.sapiens**(Bioconductor Core Team 2016a) for human) 
2. the **biomaRt** package (Durinck et al. 2005; Durinck et al. 2009) which interfaces the Ensembl genome databases in order to perform gene annotation.

这类信息包括gene symbols, gene names, chromosome names和locations, Entrez gene ID, RefSeq gene IDs，Ensembl gene IDs等.**BiomaRt** 主要用于Ensembl gene IDs,Mus.musculus包的信息来源较多，允许用户选取不同gene ID作为关键字

我们数据集使用的是Entrez gene ID，可以用Mus.musculus提取相关的gene symbols和chrmosome information

```r
gene_id <- rownames(x)
genes <- select(Mus.musculus,keys=gene_id,columns=c("SYMBOL","TXCHROM"),keytype="ENTREZID")
```

和其他的基因ID一样，Entrez gene ID也存在一对多的现象。**我们需要检查哪些gene ID是重复的，并且在处理重复前先理解source of duplication**.我们的基因注释包含28个map到多染色体的基因，（例如Gm1987就和chr4,_chr4_JH584294_random_关联，microRNA Mir5098 is associated with _chr2_, _chr5_, _chr8_, _chr11_ and _chr17_)
有两种解决方案，一个将所有重复的染色体信息的多比对gene整合到一起，或者挑选其中一个。后面的方法最方便

```r
genes <- genes[!duplicated(genes$ENTREZID),]
x$genes <- genes
```

# 数据预处理

## Transformations from the raw-scale

对于差异表达和其相关分析而言，由于文库测序深度越高会导致更高的counts，所以基因表达不会在raw counts水平进行。反之，将raw counts转换成与文库大小相关的规模(scale)是常用措施。常用的转换方法有：counts per million(CPM), log2-counts per million(log-CPM), reads per kilobase of transcript per million(RPKM) and fragments per kilobase of transcipt per million(FPKM).

在我们的分析中，CPM和log-CMP转换比较常用，尽管他们没有像RPKM和FPKM一样解释feature 长度差异。在RPKM和FPKM也能用的情况下，CPM和log-CPM值可以只用count matrix进行计算，而且对于目标比较区域特别高效。**假设不同条件下不存在isoform差异**，差异表达分析观察不同条件下的表达变化，而不是比较多个基因表达或者从绝对表达量得出结论。换句话说，基因长度对比较部分而言是固定的，任何观察上的差异都是条件改变的请过而不是由于基因长度。

此处,raw count使用**edgeR**的`cpm`转换成CMP和log-CPM值，并且log-transformation会使0.25作为预设值避免log of zero.RPKM计算也非常容易，在已知基因长度时使用`rpkm`函数。

## 移除低表达基因

所有的数据集都会包括表达和未表达的基因。虽然我们想检查哪些基因在一种条件表达而在另一个条件下不表达，但是有一些基因在所有的样本中都不表达。实际上这个数据集的19%的基因是zero counts

```r
table(rowSums(x$counts==0)==9)
```

任意条件，有生物学意义水平中都不表达的基因需要被丢弃，从而使基因的子集只含目标部分，也降低了研究差异表达下游检测的数目。
检查log-CPM值时可以发现样本中的大部分基因都没有表达或者低表达。基因超过CPM=1或log-CPM=0这个阈值时被认为是表达。基因必须在至少一个组中（或在整个实验中有三个样本）表达，用于下游分析。

尽管任何有意义的值都可以作为阈值，这里选CPM=1是因为能够很好的用于大多数数据集。这里,CPM=1说明基因表达，因为它在最低测序深度 (JMS9-P8c, library size approx. 20 million)里有最少20个counts或者在最大深度(JMS8-3, library size approx. 76 million).有76个counts.如果测序读段根据exons而不是基因或实验测序深度较低，那么CPM的阈值可以选择更低。

```r
keep.exprs <- rowSums(cpm > 1) >=3
x <-x[keep.exprs,,keep.lib.sizes=FALSE]
```

在这个标准下，大约一般的基因被移除了，用于比较前后差异的代码如下

```r
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
 den <- density(lcpm[,i])
 lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
   den <- density(lcpm[,i])
   lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
```

## 标准化基因表达分布

在样本准备或测序时，一些不是研究对象的外部因素会影响到单个样本的表达。例如第一批次处理的样品与第二批次的相比可能会有更高的表达量。这是基于所有样本应该有类似表达量范围和分布。这就是需要进行标准化。

任意描绘每个样本表达部分的作图，如密度或和箱图，在判断样本之间是否相似时都特别有用。log-CPM分布在这个数据集中的所有样本都非常相似。

**edgeR**的`calcNormFactors`使用trimmed mean of M-values(TMM)进行标准化。这里计算的标准化因子用作文库大小的刻度因子(scaling factor)。当使用DEGList-objects时，这些标准化因子会自动存放在`x$samples$norm.factors`

```r
x <- calcNormFactors(x,method="TMM")
x$samples$norm.facotrs
```

为了更好的理解标准化的影响，我们对处理后数据做一些修改，然后进行可视化

```r
x2 <- x

x2$samples$norm.factors <-1
x2$counts[,1] <- ceiling(x2#counts[,1]*0.05)
x2#counts[,2] <- x2$counts[,2]*5

par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors

## [1] 0.0547 6.1306 1.2293 1.1705 1.2149 1.0562 1.1459 1.2613 1.1170

lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")
```

## 无监督样本聚类

在我们的观念中，检验基因表达分析的一个重要探索性作图就是multi-dimensional scaling(MDS) plot或类似。该作图以无监督的方式展现样本件的相似和差异性，那么我们在进行常规测验能大致知道哪些是差异表达。理想情况下，样本会在主要的感兴趣条件下有良好的聚类，任何在组别外游离样本会被识别，伴随着sources of error或额外变异。并且技术重复相互间很近。

可以使用**limma**的`plotMDS`做这个图。第一维代表leading-fold-change，分割效果最好并且能解释数据中大部分的变异，随后维度有较少的效用，与前者垂直。当实验设计和多因素有关时，建议在不同维度中检查每一个因素。如果样本根据给定因子在任意维度都能进行聚类，这意味这个因子对表达差异有贡献，应该包括在线性模型中。另一方面，那些几乎没有影响的因子应在下游分析时排除。

在这个数据集中，样本根据实验组能良好聚类（维度1和2），然后根据泳道（维度3）。记住第一维能解释数据变异的最大部分，不同维度下值的范围会越来越小，当我们移向更高维度时。

```r
# log-CPM转换数据
lcpm <- cpm(x, log=TRUE)
# 画图设定
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")
```

**Glimma**提供了交互式MDS作图，用于探索多维

```r
glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), 
          groups=x$samples[,c(2,5)], launch=FALSE)
```

# 差异表达分析

## 构建design matrix和contrasts

本次研究的目标是了解哪些基因在三个细胞群体中表达量不同。在研究中，线性模型与数据拟合，**假设背后的数据分布正常**。
在开始前，需要设置design matrix，包含cell population和sequencing lane信息。

```r
design <- model.matrix(~0+group+lane)
colnames(desgin) <- gsub("group","",colnames(design))
```

对于一个给定的试验，通常有一些等同的方法去设置恰当的design matrix。例如`~0+group+lane`从第一个因子`group`移除 intercept,但是第二个因子`lane`依然有intercept. Understanding how to interpret the coefficients estimated in a given model is key here.我们选择第一个模型用于分析，因为设置model contrasts在`group`没有intercept时更加简明。不同细胞群体的逐个比较contrast在**limma**中用`makeContrasts`函数设置

> The intercept column contains the same value,1, for all the samples, and simply indicates that the linear model which will be built for each gene will have an intercept term.
> Intercept对应平均表达水平，当所有的试验因素都处于标准状态。

```r
contr.matrix <- makeContrasts(
    BasalvsLP = Basal-LP,
    BasalvsML = Basal-ML,
    LPvsML=LP-ML,
    levels=colnames(design))
contr.matrix
```

**limma**线性建模方法的优势在于随意调节试验的复杂度。another possibility in **limma** is to estimate correlations using `duplicateCorrelation` by specifying a `block` argument for both this function and in the `lmFit` linear modelling step.

## 移除count data的异方差性(heteroscedascity)

有文献证明，RNA-Seq count数据的方差与平均值不独立，对于raw counts和变换后log-CPM数据都是正确的。model counts使用负二项分布方法**假设**存在2项平均值-方差关系。在**limma**中，使用log-CPM值进行的线性建模被认为是正态分布，并且平均值-方差关系通过`voom`函数精确计算权重进行调节。

在DGEList-object上使用时，`voom`通过自动从x本身提取文库大小和标准化因子，将raw counts转换成log-CPM值.log-CPM值的附加标准化可以用`voom`的`normalize.method`声明

一般而言，_voom plot_表现了均值和方差间减少的趋势，由测序实验中技术差异和不同群体的样本重复的生物学差异共同引起。试验有更高的生物学差异通常有更扁平的趋势，因为方差值在高表达值中变得平滑。Experiments with low biological variation tend to result in sharp decreasing trends.

此外，_voom plot_对**上游过滤水平的可视化检查**。如果低表达基因过滤不明显，在expression scale的low end可以发现方差水平的drop，这是由samll counts引起。当这个现象发生时，我们需要返回到过滤步骤，提高阈值。

虽然样本水平的差异来源于前期的MDS作图观察，但是`voomWithQualityWeights`可以用来同时将样本水平权重和`voom`估计丰度独立权重整合一起。以Liu _al al.(2016)_为例。

```r
par(mfrow=c(1,2))
v < voom(x,design,plot=TRUE)
v

vfit <- lmFit(v,design)
vfit <- contrasts.fit(vfit,contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit,main="Final model: mean-variance trend")
```

Note that the other data frames stored within the DGEList-object that contain gene- and sample-level information, are retained in the EList-object `v` created by `voom`. The `v$genes` data frame is equivalent to `x$genes`, `v$targets` is equivalent to `x$samples`, and the expression values stored in `v$E` is analogous to `x$counts`, albeit on a transformed scale. In addition to this, the `voom` EList-object has a matrix of precision weights `v$weights` and stores the design matrix in `v$design`.

## Fitting linear models for comparisons of interest

**limma**的线性建模主要使用`lmFit`和`contrasts.fit`，原本适用于微矩阵。这两个函数可以同时用microarray和RNA-Seq数据，对每个基因拟合一个独立模型到表达值中。随后根据所有基因间的信息，进行经验贝叶斯修正，获得逐个基因变动性的精确估计。以模型残差方差对平均表达量作图。不难发现均值和方差不再具有相关性。

## 检查差异表达基因数目

为了快速查看差异表达水平，显著性上升和下调的基因可以在总结在一个表格。显著性定义为默认阈值为5%修正p值

```r
summary(decideTests(efit))
#    BasalvsLP BasalvsML LPvsML
# -1      4127      4338   2895
# 0       5740      5655   8825
# 1       4298      4172   2445
# 不难看出LP和Basal相比，4127下调，4298上调，一共8425个基因有差异
```

有些研究不只需要adjusted p-value cut-off。对于严格定义的显著性，可能还需要log-fold-changes（log-FCs)

```r
tfit <- treat(vfit,lfc=1)
dt <- decideTests(tfit)
summary(dt)
##    BasalvsLP BasalvsML LPvsML
## -1      1417      1512    203
## 0      11030     10895  13780
## 1       1718      1758    182
# 明显差异表达的基因少了
```

在多重比较下差异表达的基因可以从`decideTests`结果中提取，0：无差异，1：上调，-1：下调.`write.fit`函数可以把三个比较中的结果写到一个文件中

```r
de.common <- which(dt[,1]!=0 & dt[,2]!=0
length(de.common)
head(tfit$genes$SYMBOL[de.common],n=29)
write.fit(tfit,dt,file="results.txt")
```

使用韦恩图展示LP和ML共同改变的基因。

```r
vennDiagram(dt[,1:2],circle.col=c("turquoise","salmon"))
```

## 从头到底检查独立DE基因。

高DE基因可以使用`topTreat`对于使用`treat`的结果而言。（或者是`topTable`对`eBayes`结果）。默认`topTreat`根据adjusted p-values从低到高安排基因和相关信息,log-FC, average log-CPM, moderated t-statistic, raw and adjusted p-values。

```r
basal.vs.lp <- topTreat(tfit,coef=1,n=Inf)
basal.vs.ml <- topTreat(tfit,coef=2,n=Inf)
head(basal.vs.lp)
head(basal.vs.ml)
```

## 差异表达结果可视化

为了可视化总结所有基因结果，mean-differenct plots根据线性模型拟合相对平均log-CPM值显示了线性log-FC。该图可以使用`plotMD`产生，高亮那些差异表达基因。

```r
plotMD(tfit,column=1,status=dt[,1],main=colnames(tfit)[1],xlim=c(-8,13))
```

**Glimma**的`glMDPlot`提供交互式均值-差异作图扩展上个函数。

```r
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         id.column="ENTREZID", counts=x$counts, groups=group, launch=FALSE)
```

上面的图包括咋任何条件下都表达的基因，如维恩图，，或者单独观察基因
**Heatmap**允许我们观察部分基因。提供一种好用的角度，在兼顾总体研究的同事，关注了单个基因。在关注上千个基因的同时也保证了解析度。

**gplots**包的`heatmap.2`根据basal vs LP创建了前100DE基因。热图正确的根据细胞类型聚类了样本，重新排布了基因顺序组成了类似表达模块。