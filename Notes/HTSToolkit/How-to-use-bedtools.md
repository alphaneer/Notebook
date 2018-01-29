---
title: 如何使用bedtools处理Ranges数据
tags: 工具笔记, HTS
notebook: 工具笔记
---

# 如何使用bedtools处理 _Rang_ 数据

## 什么是 _Range_ 数据

参考基因组表示的是一种坐标系统，比如说某一个物种基因组大小为100bp，那么他参考基因组就可以表示为[1,100], 之后就可以用任意[x,y]表示这条参考基因组上的位置，这就是一种范围信息，X-Y这段区域可能是外显子，也可能是内含子，可能是编码区，也可能是基因间区，也有可能是一个测序结果。

因此 _Range_ 数据是生信数据比较常见的存放形式，比如说BED/BAM/BCF/和GFF/BFF/SAM/VCF/，前者以0为始，后者以1为始。

为了操作这种 _Range_ 数据，Bioconductor在R语言中定义了两个重要的对象，IRange和GenomicRanges，后者仅存放'start','end','width'是后者的基础。后者才能真正存放基因组 _Range_ 数据。

这一篇不介绍如何在R语言操作 _Range_ 数据，而是介绍bedtools这款号称基因组 _Range_ 数据分析的瑞士军当，当时的口号是一款取代10个生信分析师的工具。

Bedtools能够对基因组 _Range_ 数据进行**交**，**并**，**补**，**计数**等简单操作，也能和Unix命令行结合起来完成更加复杂的任务。

## bed格式简介

在正式介绍bedtools之前，需要先介绍一下BED格式。根据USCSC基因组浏览器的描述，BED格式能够非常简洁的表示基因组特征和注释，尽管BED格式描述中定义了12列，但是仅仅只有3列必须，因此BED格式按照列数继续细分为BED3,BED4,BED5,BED6,BED12。

BED12定义的12列分别为：chrom, start, end, name(BED代表的特征名),score(范围为0~1000，可以是pvalue, 或者是字符串,如"up"), strand(正负链), thickstart, thickednd(额外着色位置, 比如说表示外显子), itemRgb(RGB颜色,如255,0,0), blockCount(区块数量, 如外显子), blockSizes(由逗号隔开的区块大小), blockStarts(由逗号隔开的区块起始位点)。

知道了BED12后，就可以对BED的细分格式进行举例说明

- BED3：`chr1          11873   14409`
- BED4: `chr1  11873  14409  uc001aaa.3`
- BED5: `chr1 11873 14409 uc001aaa.3 0`
- BED6: `chr1 11873 14409 uc001aaa.3 0 +`
- BED12: `chr1 11873 14409 uc001aaa.3 0 + 11873 12000 123,123,123 3 354,109,1189, 0,739,1347,`

![BED12效果](http://oex750gzt.bkt.clouddn.com/18-1-24/37066020.jpg)

除了官方的BED定义外，BEDtools定义了BEDPE用来表示基因组不连续的特征，比如说结构变异或者双端测序的reads。在定义中10列是必须的，为chrom11, start1, end1, chrom2, start2, end2, name, score, strand1, strand2。 这之后可以增加任意多的其他列。

其他BEDtools支持的格式说明：

- GFF3: <https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md>
- SAM/BAM: <http://samtools.github.io/hts-specs/SAMv1.pdf>
- VCF: <http://samtools.github.io/hts-specs/VCFv4.2.pdf>

## Bedtools工具介绍

bedtools的功能非常强大，试图解决你所遇到的所有和基因组位置运算的问题以及周边问题：基因组运算，多文件比较， PE数据操作，格式转换， Fasta数据操作， BAM工具， 统计学相关工具，其他小工具

其中最重要的选项是`--help`，一个强大的工具提供了许多参数，需要勤读帮助文档。bedtools的官方文档写的非常优秀，绝大部分工具都以图解的方式形象的说明了每个参数的可能效果。因此我写这篇文章的目的就是让迫使自己去熟悉所有的工具而已。

### 核心：基因组运算

所谓的基因组运算，就是看看看自己手头拿到的区域和你感兴趣的区域的关系如何。bedtools提供了如下工具做一系列你想到或者你想不到的事情。

集合运算：

- intersect: **交集**，也就是看两个区域的重叠
- window: 和intersect类似也是求交集，但是还会看上下游一定区间内是否有重叠。 比如说看某个基因的上下游是否有peak
- complement: **补集**，根据一个已有区域得到另一个不重叠的区域
- substract: **差集**， 求一个区域减去另一个区域后的结果
- merge: **合集**，将相邻的区间合并成一个
- cluster: 类似于merge，但是不做合并

区间统计：

- closest: 找到和目标最近的特征区域。比如说ChIP-seq得到的peaks可以根据“距离其最近的基因就是他调控的基因”的假设进行注释。
- coverage: 计算某个区间的覆盖情况
- genomecov: 计算全部基因组的覆盖情况
- map: 在某个区间的应用其他函数， 比如说求和(sum), 均值(mean)
- annotate: 从其他一系列文件中统计给定区间的覆盖情况

区间工具：

- flank: 根据已有特征区间，得到两翼位置的新区间
- slop: 根据已有特征区间，向外衍生
- shift: 特征区域整体偏移一定位置
- shuffle: 从基因组区间中随机选择某些位置
- sample: 从已有的区域随机挑选
- makewindows: 从基因组上创建等长区间

### 其他

bedtools的核心工具就是上面几个，剩下的都算是小轮子，解决了你手上轮子不够多的烦恼。

Fasta相关：

- getfasta: 根据提供的区间提取序列
- maskfasta: 根据提供的区间对序列进行遮盖
- nuc: 统计某个区间的碱基含量

BAM格式相关：

- 格式转换：bamtobed, bedtobam, bamtofastq,bedpetobam,bed12tobed6
- multicov: 计算某些位点上不同BAM文件覆盖
- tag: 计算一个BAM在不同位置上覆盖

PE文件操作：

- pairtobed: 专门用于看其他格式在BEDPE格式上重叠情况
- pairtopair: 比较两个BEDPE格式的重叠情况

统计学工具：主要是用以不同的统计学方法来衡量两个区间的相似度，有三种： jaccard, reldist, fisher

除了以上，还有一些更加有趣的小工具，比如说`igv`可以创建IGV自动截图的运行脚本，`links`可以构建能在UCSC基因组上打开的链接等。

## bedtools使用案例【待续】
