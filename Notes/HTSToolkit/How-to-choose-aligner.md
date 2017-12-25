---
title: 如何比对高通量测序结果
tags: 高通量
notebook: 工具笔记
---
# 高通量数据分析初步：序列比对

在过去的十几年里，随着高通量测序(HTS)成本降低，出现了各种测序概念, DNA-Seq, ChIP-Seq, RNA-Seq, BS-Seq覆盖了研究领域的方方面面。随着而来的问题是，如何把这些短片段**快速且准确**地回贴到参考基因组上。

## 为什么要需要高通量测序比对工具

## align vs map

## 如何挑选合适的短读比对工具

2012年 _Bioinformatics_ 有一篇文章^[Tools for mapping high-throughput sequencing data ]综述了目前高通量数据的比对软件，并且建立主页<https://www.ebi.ac.uk/~nf/hts_mappers/>罗列并追踪目前的比对软件。

![](http://oex750gzt.bkt.clouddn.com/17-12-23/79073456.jpg)

尽管看起来有那么多软件，但是实际使用就那么几种，BWA(傲视群雄), TopHat(尽管官方都建议用HISAT2，还是那么坚挺), SOAP(架不住华大业务多)。 Blat和Mummer3只能短读长序列比对功能，应该不能归为高通量一行。i

## 两个工具：BWA 和 Bowtie2

## 如何比较工具的好坏
