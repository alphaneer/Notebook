---
title: 纯二代测序从头组装动植物基因组
tags: 组装, Bioinformatics
notebook: 软件工具
---
# 基因组组装

基因组组装一般分为三个层次，contig, scaffold和chromosomes. contig表示从大规模测序得到的短读(reads)中找到的一致性序列。组装的第一步就是从短片段(pair-end)文库中组装出contig。进一步基于不同长度的大片段(mate-pair)文库，将原本孤立的contig按序前后连接，其中会调整contig方向以及contig可能会存在开口(gap,用N表示)，这一步会得到scaffolds,就相当于supercontigs和meatacontigs。最后基于遗传图谱或光学图谱将scaffold合并调整，形成染色体级别的组装(chromosome).

目前基于二代测序的组装存在挑战：

- 全基因组测序得到的短读远远小于原来的分子长度
- 高通量测序得到海量数据会增加组装的计算复杂性，消耗更高的计算资源
- 测序错误会导致组装错误，会明显影响contig的长度
- 短读难以区分基因组的重复序列
- 测序覆盖度不均一，会影响统计检验和结果结果诊断

上述的问题可以尝试从如下角度进行解决

- 短读长度：可以通过提供更多样本，并且建库时保证位置足够随机
- 数据集大小: 使用K-mers算法对数据进行组装。assembler不再搜寻overlap，而是搜索具有相同k-mers的reads。但是k-mer算法相比较overlap-based算法，灵敏度有所欠缺，容易丢失一些true overlaps。关键在于定义K。**注**: K-mer表示一条序列中长度为k的连续子序列,如ABC的2-mer为AB,BC
- 测序错误: 必须保证测序结果足够正确, 如提高质量控制的标准
- 基因组重复区： 测序深度要高，结果要正确。如果repeat短于read长度，只要保证有足够多且特异的read。如果repeat长于read，就需要paired ends or “mate-pairs”
- 覆盖度不均一： 提高深度，保证随机
- 组装结果比较：contig N50, scaffold N50, BUSCO

![什么叫做N50](http://bmpvieira.github.io/assembly14/img/n50.png)

## 二代数据从头组装的主流工具

基因组组装的组装工具主要分为三类：基于贪婪算法的拼接方法，基于读序之间的重叠序列(overlapped sequence)进行拼接的OLC(Overlap-Layout-Consensus)拼接方法和基于德布鲁因图的方法，这三种方法或多或少基于图论。第一种是最早期的方法，目前已被淘汰，第二种适用于**一代测序**产生长片段序列，第三种是目前二代测序组装基因组的工具的核心基础。

![几种比较复杂的图](http://oex750gzt.bkt.clouddn.com/18-3-1/92185065.jpg)

根据"GAGE: A critical evaluation of genome assemblies and assembly algorithms"以及自己的经验，目前二代数据比较常用的工具有Velvet, ABySS, AllPaths/AllPaths-LG, Discovar, SOAPdenovo, Minia, spades。这些工具里，ALLPaths-LG是公认比较优秀的组装工具，但消耗内存大，并且要提供至少两个不同大小文库的数据, SOAPdenovo是目前使用率最高的工具(华大组装了大量的动植物基因组)。工具之间的差别并没有想象的那么大，在物种A表现一般的工具可能在物种B里就非常好用，因此要多用几个工具，选择其中最好的结果。

## 组装的基本流程

在正式组装之前，需要先根据50X左右的illumina测序结果对基因组进行评估，了解基因组的大小，重复序列含量和复杂度。基于这些信息，确定后续策略以及是否真的需要对该物种进行测序。

当你拿到测序数据后，就可以按照如下几步处理数据。第一步是**数据质控控制**，这一步对于组装而言非常重要，处理前和处理后的组装结果可能会天差地别；第二步，根据经验确定**起始参数**，如K-mer和覆盖率；第三步，使用不同软件进行组装；第四步，评估组装结果，如contig N50, scaffold N50, 判断是否需要修改参数重新组装。

## 数据准备

这里使用来自于[GAGE](http://gage.cbcb.umd.edu/data/index.html)的金黄色葡萄球菌 _Staphylococcus aureusa_ 数据进行练习。一方面数据量小，服务器能承受并且跑得快，另一方面本身基因组就组装的不错，等于是考完试能够自己对答案。

```bash
mkdir Staphylococcus_aureus && cd Staphylococcus_aureus
mkdir genome
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz > genome/Saureus.fna.gz
mkdir -p raw-data/{lib1,lib2}
curl http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/frag_1.fastq.gz > raw-data/lib1/frag_1.fastq.gz
curl http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/frag_2.fastq.gz > raw-data/lib2/frag_2.fastq.gz
curl http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/shortjump_1.fastq.gz > raw-data/lib2/shortjump_1.fastq.gz
curl http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/shortjump_1.fastq.gz > raw-data/lib2/shortjump_1.fastq.gz
```

**双链特性**：forward sequence of a read may overlap either the forward sequence or the reverse complement sequence of other reads
**回文序列**：k-mers长度为奇数
**测序错误**：  数据预处理， graph edges which represent  a higher number of K-mers are more highly trusted, sequence alignment algorithms are used to collapse highly similar paths
**重复**：inverted repeats, tandem repeats, inexact repeats, and nested repeats， 增加k-mer长度，或mate-pair 测序
