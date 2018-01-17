---
title: Bioconductor数据类型简介
tags: R, bioconductor
notebook: 工具笔记
---

# Introduction to Bioconductor for Sequence Data

在院士组轮转的时候，由于没人安排我去做什么，我也不懂怎么和别人交流，于是那些日子就在翻译Biocondutor的教程中度过。本文为我学习Bioconductor所写的第一篇学习笔记。说是学习笔记，差不多等于全文翻译了bioconductor的 **sequence Analysis**一章。基本上高通量数据分析到了最后都会用到Bioconductor的包，甚至是前期都可能会用到，

## 什么是Bioconductor

官方的简介是：
> Bioconductor provides tools for the analysis and comprehension of high-throughput genomic data. Bioconductor uses the R statistical programming language, and is open source and open development. It has two releases each year, [1296 software packages](http://www.bioconductor.org/packages/release/BiocViews.html#___Software), and an active user community. Bioconductor is also available as an [AMI](http://www.bioconductor.org/help/bioconductor-cloud-ami/) (Amazon Machine Image) and a series of [Docker](http://www.bioconductor.org/help/docker/) images.

简单的说Bioconductor就是以R语言为平台的一个高通量基因组数据的分析工具。那么下面就算正文环节

## 简介

Biconductor为高通量基因组数据提供了分析和理解方法。我们拥有大量包可对大量数据进行严格的统计分析，while keeping techological aritiacts in mind. Bioconductor能进行多种类型分析：

- 测序： RNASeq, ChIPSeq, variants, copy number...
- 微矩阵： 表达，SNP
- Domain specific analysis: Flow cytometry, Proteomics...

对于这些分析，可用简单的导入并使用多种序列相关文件类型，包括,fasta, fastq,BAM,gtf,bed,和wig文件等等.Bionconductor支持导入，常见和高级的序列操作，triming,transformation, and alignment(质量评估)

## Sequencing Resources

不同阶段，不同数据格式会用到的包

![](http://oex750gzt.bkt.clouddn.com/18-1-16/62012287.jpg)

- **IRange, GenomicRanges**：基于range的计算，数据操作，常见数据表示。**Biostring**: DNA和氨基酸序列表示，比对和模式匹配和大规模生物学序列的数据操作。**ShortRead**处理FASTQ文件
- **Rsamtools, GenomicAlignments**用于比对read(BAM文件)的I/O和数据操作 **rtracklayer**用于导入和导出多种数据格式(BED,WIG,bigWig,GTF,GFF)和UCSC genome browser tracks操作
- **BSgenomes**： 用于获取和操作规划的全基因组。**GenomicFeatures**常见基因组之间序列特征注释，**biomaRt** 获取Biomart数据库
- **SRAdb**用于从Sequnce Read Archive(SRA)查找和获取数据。

Bioconductor包通过biocViews进行管理，有些词条位于Sequencing，其他则在其他条目和代表性的包下，包括

* [RNASeq](http://bioconductor.org/packages/biocViews.html#__RNASeq), e.g., _[edgeR](http://bioconductor.org/packages/edgeR)_, _[DESeq2](http://bioconductor.org/packages/DESeq2)_, _[edgeR](http://bioconductor.org/packages/edgeR)_, _[derfinder](http://bioconductor.org/packages/derfinder)_, and _[QuasR](http://bioconductor.org/packages/QuasR)_.

* [ChIPSeq](http://bioconductor.org/packages/biocViews.html#__ChIPSeq), e.g.,_[DiffBind](http://bioconductor.org/packages/DiffBind)_, _[csaw](http://bioconductor.org/packages/csaw)_, _[ChIPseeker](http://bioconductor.org/packages/ChIPseeker)_, _[ChIPQC](http://bioconductor.org/packages/ChIPQC)_.

* [SNPs](http://bioconductor.org/packages/biocViews.html#__SNP) and other variants, e.g., _[VariantAnnotation](http://bioconductor.org/packages/VariantAnnotation)_, _[VariantFiltering](http://bioconductor.org/packages/VariantFiltering)_, _[h5vc](http://bioconductor.org/packages/h5vc)_.

* [CopyNumberVariation](http://bioconductor.org/packages/biocViews.html#__CopyNumberVariation) e.g., _[DNAcopy](http://bioconductor.org/packages/DNAcopy)_, _[crlmm](http://bioconductor.org/packages/crlmm)_, _[fastseg](http://bioconductor.org/packages/fastseg)_.

* [Microbiome](http://bioconductor.org/packages/biocViews.html#__Microbiome) and metagenome sequencing, e.g., _[metagenomeSeq](http://bioconductor.org/packages/metagenomeSeq)_, _[phyloseq](http://bioconductor.org/packages/phyloseq)_, _[DirichletMultinomial](http://bioconductor.org/packages/DirichletMultinomial)_.

## Ranges Infrastructure

许多其他的Bioconductor包高度依赖于IRanges/GenomicRanges所提供的低层数据结构

```R
library(GenomicRanges)
GRanges(seqnames=Rle(c('chr1', 'chr2', 'chr3'), c(3, 3, 4)),
      IRanges(1:10, width=5), strand='-',
      score=101:110, GC = runif(10))
```

GenomicRanges使得我们可以将染色体坐标的一个范围和一段序列名（如chromosome)以及正负链关联起来。这对描述数据(aligned reads坐标,called ChIP peaks 或copy number variants)和注释(gene models, Roadmap Epigenomics regularoty elements, dbSNP的已知临床相关变异)

GRanes对象，用来储存基因组位置和关联注释的向量。
Genomic ranges基本都是通过导入数据 (e.g., via`GenomicAlignments::readGAlignments()`)或注释(e.g., via `GenomicFeatures::select()` or `rtracklayer::import()` of BED, WIG, GTF, and other common file formats).进行创建

```R
help(package="GenomicRanges")
vignette(package="GenomicRanges")
```

GRanges常用操作包括findOverlaps， nearest等

## FASTA文件中的DNA/氨基酸序列

**Biostrings**类用于表示DNA或氨基酸序列

```R
library(Biostrings)
d <- DNAString("TTGAAAA-CTC-N")
length(d)
```

使用AnnotationHub在Ensembl上下载 _Homo sapiens_ cDNA序列，文件名为‘Homo_sapiens.GRCh38.cdna.all.fa’

```R
library(AnnotationHub)
ah <- AnnotationHub()
ah2 <- query(ah, c("fasta","homo sapiens","Ensembl"))
fa <- ah2[["AH18522]]
fa
```

使用Rsamtools打开fasta文件读取里面的序列和宽度记录

```R
library(Rsamtools)
idx <- scanFaIndex(fa)
idx

long <- idx[width(idx)>82000]
getSeq(fa,param=long)
```

**BSgenome**提供全基因组序列，

```R
library(BSgenome.Hsapiens.UCSC.hg19)
chr14_range <- GRanges("chr14",IRanges(1,seqlengths(HSapiens)["chr14"]))
chr12_dna <- getSeq(Hsapiens,chr14_range)
letterFrequency(chr14_dna,‘GC",as.prob=TRUE)
```

## FASTQ文件中的Reads

**ShortRead**用于fastq文件.**BiocParallel**用于加速任务进程

```R
##1. attach ShortRead and BiocParallel
library(ShortRead)
library(BiocParallel)
## 2. create a vectior of file paths
fls <- dir("~/fastq",pattern="*fastq",full=TRUE)
## 3. collect statistics
stats0 <- qa(fls)
## 4. generate and browser the report
if (interactive())
    browseURL(report(stats))
```

## BAM中的Aligned Reads

**GenomicAlignments**用于输入比对到参考基因组的reads
下一个案例中，我们会读入一个BAM文件,and specifically read in reads / supporting an apparent exon splice junction/ spanning position 19653773 of chromosome 14.

 _[RNAseqData.HNRNPC.bam.chr14_BAMFILES](http://bioconductor.org/packages/RNAseqData.HNRNPC.bam.chr14_BAMFILES)_ 内有8个BAM文件。我们只是用第一个BAM文件。我们会载入软件包和数据包，为目标区域构建一个GRanges对象，然后使用 _summarizeJunctions()_ 寻找目标区域的read

```R
# 1. 载入软件包
library(GenomicRanges)
library(GenomicAlignments)
# 2. 载入样本数据
library('RNAseqDta.HNRNPC.bam.chr14')
bf <- BamFile(RNAseqData.HNRNPC.bam.chr14_BAMFILES[[1),asMates=TRUE)
# 3. 定义目标区域
roi <- GRanges("chr14", IRanges(19653773,width=1)）
# 4. alignments, junctons, overlapping our roi
paln <- readAlignmentsList(bf)
j <- summarizeJunctions(paln,with.revmap=TRUE)
j_overlap <- j[j %over% roi]

# 5.supporting reads
paln[j_overlap$revmap[[1]]]
```

## Called Variants from VCF files

VCF(Variant Call Files)描述了SNP和其他变异。该文件包括元信息行，含有列名的header line和众多的数据行，每一行都包括基因组上位置信息和每个位置上可选基因型信息。
数据通过 _VariantAnnotation_ 的readVcf()读入

```R
library(VariantAnnotation)
fl <- system.file("extdata","chr22.vcf.gz",package=“VariantAnnotation")
vcf <- readVcf(fl,"hg19")
```

## Genome Annotations from BED, WIG, GTF etc files

[rtracklayer](http://bioconductor.org/packages/rtracklayer)能够导入和导出许多常见的文件类型，BED,WIG,GTF,还能查询和导航UCSC genome Browser
_rtracklayer_包含测序所用BED文件

```R
library(rtracklayer)
test_path <- system.file("test",package = "rtracklayer")
test_bed <- file.path(test_path,"test.bed")

test <- import(test_bed,format="bed")
test
```

_[AnnotationHub](http://bioconductor.org/packages/AnnotationHub)_ 同样包括多种基因组注释文件（BED,GTF,BigWig)，使用rtracklayer import()。更多有关AnnotationHub的介绍需要看 [Annotation workflow](http://bioconductor.org/help/workflows/annotation/Annotation_Resources/) and [AnnotationHub HOW TO vignette](http://bioconductor.org/packages/devel/bioc/vignettes/AnnotationHub/inst/doc/AnnotationHub-HOWTO.html)。