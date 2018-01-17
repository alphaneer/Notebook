---
title: 基因组变异位点注释
tags: R, bioconductor
notebook: 工具笔记
---
# 基因组变异位点注释

安装工作流程所需的`biconductor`包

```r
source("http://bioconductor.org/workflows.R")
workflowInstall("variants")
```

## 背景

VariantAnnotation包能够有效的从Variant Calling Format(VCF)文件读取部分或所有内容。
这些文本文件包括元信息行(meta-information lines)，标题行(header line)和数据行(data lines)，其中数据行每一行都含有基因组位置信息。这类格式同样包含每个位置上样本的基因型信息。更多该文件相关的信息可以看[VCF specs](http://samtools.github.io/hts-specs/VCFv4.2.pdf)

## 配置

本文所介绍的工作流程需要一些Biocondutor的包，下面几节会仔细介绍每个包的具体用法。

```r
library(VariantAnnotation)
library(cgdv17)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(PolyPhen.Hsapiens.dbSNP131)
```

可以用`biocLite`安装那些未安装的包

```r
source("https://bioconductor.org/biocLite.R")
biocLite("mypackage")
```

## 探索TRPV基因家族的变异位点

本工作流程着眼于17号染色体上Transient Receptor Potential Vanilloid (TRPV)基因家族的变异位点。样本数据来自于Bioconductor的cgdv17实验数据包，内部包含46个17号染色体上的完整的基因组多样性面板数据（pannel data).如果想知道这些数据是如何组织的信息，可以查看包的小品文。

```r
browseVignettes("cgdv17")
```

我们所使用的包中的VCF文件，是CEU群体其中一个17号染色体的子集。

```r
library(VariantAnnotation)
library(cgdv17)
file <- system.file("vcf", "NA06985_17.vcf.gz", package = "cgdv17")
```

### 检查VCF文件的标题数据

为了大致了解该文件有哪些数据，我们可以查看标题部分。`scanVcfHeader()`解析文件的标题部分，将解析的内容存入`VCFHeader`对象，然后就可以使用`info()`和`geno()`存取器（accessor）提取字段特定(field-specific)数据.

```r
hdr <- scanVcfHeader(file)
info(hdr)
geno(hdr)
```

由下可知，VCF中的变异比对到NCBI构建的基因组GRCh37.

```r
meta(hdr)$META
## DataFrame with 5 rows and 1 column
##                       Value
##                 <character>
## fileformat          VCFv4.1
## fileDate           20111102
## source     masterVar2VCFv40
## reference    build37.fa.bz2
## phasing             partial
```

### 将基因符号(symbol)变为基因ID

使用`org.Hs.eg.db`包将基因符号转为基因ID。

```r
library(org.Hs.eg.db)
## 定义哪些基因
genesym <- c("TRPV1", "TRPV2", "TRPV3")
## 从org.Hs.eg.db挑选这些符号的ID,select如何使用见之前文章
geneid <- select(org.Hs.eg.db, keys=genesym, keytype="SYMBOL",
                 columns="ENTREZID")
geneid
##   SYMBOL ENTREZID
## 1  TRPV1     7442
## 2  TRPV2    51393
## 3  TRPV3   162514
```

### 创建基因范围(gene ranges)

我们使用USCS的hg19已知基因轨道（hg19 known gene track)识别TRPV基因范围。这些基因范围最终会根据VCF文件提取变异位点。
载入注释包

```r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
```

我们的VCF已经比对到NCBI的基因组，并且已知基因轨道来自于UCSC。这些机构对染色体有不同的命名传统。为了在匹配（match)或者重叠(overlap)操作用到这些数据,染色体命名方式（或者叫seqlevels)需要匹配。我们会修改txdb以匹配VCF文件.

```r
txdb <- renameSeqlevels(txdb, gsub("chr", "", seqlevels(txdb)))
txdb <- keepSeqlevels(txdb, "17")
```

根据基因创建转录本列表

```r
txbygene = transcriptsBy(txdb, "gene")
```

为TRPV基因创建基因范围

```r
gnrng <- unlist(range(txbygene[geneid$ENTREZID]), use.names=FALSE)
names(gnrng) <- geneid$SYMBOL
```

### 提取变异位点子集

`ScanVcfParam`对象用于提取数据子集。该对象能够指定基因组坐标（范围）或单独的VCF元素。提取范围(vs 字段)需要一个tabix索引。使用`?indexTabix`查看细节。

```r
param <- ScanVcfParam(which = gnrng, info = "DP", geno = c("GT", "cPd"))
param
## 从VCF文件提取TRPV范围
vcf <- readVcf(file, "hg19", param)
## 查看VCF对象 'fixed', 'info' and 'geno'
vcf
head(fixed(vcf))
geno(vcf)
```

### 基因模型的变异位点位置

`locateVariants`根据基因结构（例如exon, utr, splice site等）判断变异位点的位置。我们使用之前加载的`TxDb.Hsapiens.UCSC.hg19.knownGene`包内的基因模型。

```r
## 使用'region'参数定义目标区域
## 详细内容见？locateVariants
cds <- locateVariants(vcf, txdb, CodingVariants())
five <- locateVariants(vcf, txdb, FiveUTRVariants())
splice <- locateVariants(vcf, txdb, SpliceSiteVariants())
intron <- locateVariants(vcf, txdb, IntronVariants())

all <- locateVariants(vcf, txdb, AllVariants())
```

CDS的每一行都代表一个变异位点-转录本匹配，因此一行变异位点对应多行也是可以的。如果我们对基因中心的问题感兴趣，数据就可以根据基因进行描述性分析，而不用考虑转录本。

```r
## variants是一对多关系么
table(sapply(split(mcols(all)$GENEID, mcols(all)$QUERYID), 
      function(x) length(unique(x)) > 1))

## 根据基因分析变异位点数量
idx <- sapply(split(mcols(all)$QUERYID, mcols(all)$GENEID), unique)
sapply(idx, length)


## 根据基因分析变异位点位置
sapply(names(idx),
    function(nm) {
        d <- all[mcols(all)$GENEID %in% nm, c("QUERYID", "LOCATION")]
        table(mcols(d)$LOCATION[duplicated(d) == FALSE])
    })
```

## 非同义突变位点的氨基酸编码改变

可用`predictCoding`函数得到非同义变异的氨基酸改变。`BSgenome.Hsapiens.UCSC.hg19`包用作参考等位基因的源。变异的等位基因由使用者提供。

```r
library(BSgenome.Hsapiens.UCSC.hg19)
seqlevelsStyle(vcf) <- "UCSC"
seqlevelsStyle(txdb) <- "UCSC"
aa <- predictCoding(vcf, txdb, Hsapiens)
```

`predictCoding`仅仅返回编码变异位点的结果。与`locateVariants`一样，每个变异位点-转录匹配项的输出都有一行，因此每个变异位点可以有多行。

```r
## 是否有多个突变位点位于同一基因
table(sapply(split(mcols(aa)$GENEID, mcols(aa)$QUERYID), 
        function(x) length(unique(x)) > 1))

## FALSE
##    17

## 根据基因总结变异位点数量
idx <- sapply(split(mcols(aa)$QUERYID, mcols(aa)$GENEID, drop=TRUE), unique)
sapply(idx, length)

## 162514  51393   7442
##      6      3      8

## 根据基因总结变异位点后果
sapply(names(idx),
       function(nm) {
           d <- aa[mcols(aa)$GENEID %in% nm, c("QUERYID","CONSEQUENCE")]
           table(mcols(d)$CONSEQUENCE[duplicated(d) == FALSE])
       })

##                162514 51393 7442
## nonsynonymous       2     0    2
## not translated      1     0    5
## synonymous          3     3    1
```

当`predictCoding`调用时,变异位点“not translated”在抛出的警告进行说明。 在varAllele中缺少varAllele或“N”的变异位点不会被翻译。 如果varAllele替换已经导致了移位，则后果将是“frameshift”。 有关详细信息，请参阅`?predictCoding`

## 使用ensemblVEP包进行注释

`ensemblVEP`包能够访问在线Ensembl Variant Effect Predictor (VEP tool)。VEP工具输出的已知或者未知变异位点的功能后果预测，通过序列本体论(Sequence Ontology)或Ensembl报告。可选输出有Regulatory region consequences, HGNC, Ensembl protein identifiers, HGVS, co-located variants。 `ensemblVEP()`接受VCF文件名，在R工作环境中返回一个磁盘上的VCF或者GRanges.
加载ensemblVEP:

```r
library(ensemblVEP)
```

ensemblVEP的`file`参数必须是硬盘上的VCF文件。我们之后会以TRPV变异位点写出VCF对象，并且上传到ensemblVEP.
以`tempfile()`写出一个VCF对象到vcf文件中

```r
dest <- tempfile()
writeVcf(vcf, dest)
```

使用只含有TRPV变异位点的文件和自定义`VERPParam`对象调用`ensemblVEP`Call ensemblVEP

```r
gr <- ensemblVEP(file = dest)
head(gr, 3)
```

来自VEP工具的数据作为元数据列返回给GRanges（'gr'对象）。 您可以通过设置`VEPParam`中的运行时选项来进一步控制返回哪些字段。

```r
VEPParam()
```

每个选项，例如“input”，“basic”，“cache”等具有可以设置的多个字段。

```r
basic(VEPParam())
```

更多运行选项的信息可以在`VEPParam`的man page或者
 [Ensembl VEP web site](http://uswest.ensembl.org/info/docs/tools/vep/script/vep_options.html)中找到.

## 了解更多信息

想要更多了解`VariantAnnotation`的使用方法，有以下方法：

```r
help(package="VariantAnnotation")
?predictCoding
browseVignettes(package="VariantAnnotation")
```
