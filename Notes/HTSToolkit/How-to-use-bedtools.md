---
title: 如何使用bedtools处理Ranges数据
tags: 工具笔记, HTS
notebook: 工具笔记
---

# 如何使用bedtools处理_Rang_数据

## 什么是_Range_数据

参考基因组表示的是一种坐标系统，比如说某一个物种基因组大小为100bp，那么他参考基因组就可以表示为[1,100], 之后就可以用任意[x,y]表示这条参考基因组上的位置，这就是一种范围信息，X-Y这段区域可能是外显子，也可能是内含子，可能是编码区，也可能是基因间区，也有可能是一个测序结果。

因此_Range_数据是生信数据比较常见的存放形式，比如说BED/BAM/BCF/和GFF/BFF/SAM/VCF/，前者以0为始，后者以1为始。

为了操作这种_Range_数据，Bioconductor在R语言中定义了两个重要的对象，IRange和GenomicRanges，后者仅存放'start','end','width'是后者的基础。后者才能真正存放基因组_Range_数据。

这一篇不介绍如何在R语言操作_Range_数据，而是介绍bedtools这款号称基因组_Range_数据分析的瑞士军当，当时的口号是一款取代10个生信分析师的工具。

Bedtools能够对基因组_Range_数据进行**交**，**并**，**补**，**计数**等简单操作，也能和Unix命令行结合起来完成更加复杂的任务。

## Bedtools使用介绍

### Bedtools模块介绍

bedtools的功能非常强大，分为如下几个模块

- 基因组运算
- 多文件比较
- PE数据操作
- 格式转换
- Fasta数据操作
- BAM工具
- 统计学相关工具
- 小工具

其中最重要的选项是`--help`，一个强大的工具提供了许多参数，需要勤读帮助文档。如果你不知道bedtools能用来干嘛，你可以看下别人的操作。

* [Coverage analysis for targeted DNA capture](http://gettinggeneticsdone.blogspot.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html). Thanks to [Stephen Turner](https://twitter.com/genetics_blog).
* [Measuring similarity of DNase hypersensitivity among many cell types](https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md#bp6--measuring-dataset-similarity)
* [Extracting promoter sequences from a genome](http://www.biostars.org/p/17162/)
* [Comparing intersections among many genome interval files](http://www.biostars.org/p/13516/)
* [RNA-seq coverage analysis](http://www.cureffi.org/2013/11/18/an-mrna-seq-pipeline-using-gsnap-samtools-cufflinks-and-bedtools/). Thanks to [Erik Minikel](https://twitter.com/cureffi).
* [Identifying targeted regions that lack coverage](https://twitter.com/aaronquinlan/status/421786507511205888). Thanks to [Brent Pedersen](https://twitter.com/brent_p).
* [Calculating GC content for CCDS exons](http://www.biostars.org/p/47047/).
* [Making a master table of ChromHMM tracks for multiple cell types](https://gist.github.com/arq5x/3138599).

### 基因组运算

后续所用到的数据可按照如下方法下载

```bash
curl -O https://s3.amazonaws.com/bedtools-tutorials/web/maurano.dnaseI.tgz
curl -O https://s3.amazonaws.com/bedtools-tutorials/web/cpg.bed
curl -O https://s3.amazonaws.com/bedtools-tutorials/web/exons.bed
curl -O https://s3.amazonaws.com/bedtools-tutorials/web/gwas.bed
curl -O https://s3.amazonaws.com/bedtools-tutorials/web/genome.txt
curl -O https://s3.amazonaws.com/bedtools-tutorials/web/hesc.chromHmm.bed
tar xf dnaseI.tgz
```

#### 求交集

`intersect`（交集）比较多个BED/BAM/VCF/GFF文件，以多种方式计算他们的交集。

![intersect](http://bedtools.readthedocs.org/en/latest/_images/intersect-glyph.png)

![merge](http://bedtools.readthedocs.org/en/latest/_images/merge-glyph.png)

![complement](http://bedtools.readthedocs.org/en/latest/_images/complement-glyph.png)

![genomecov](http://bedtools.readthedocs.org/en/latest/_images/genomecov-glyph.png)