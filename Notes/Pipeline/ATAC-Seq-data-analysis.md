---
title: ATAC-Seq数据分析
author: xuzhougeng
tag: ATAC-Seq, epigenetic
notebook: 分析流程
---
# ATAC-Seq 数据分析

## 课程要求

1. 目前生信培训很多，但是还没有ATAC-seq相关的教程，如果有，你这是第一篇，以后你就是个有作品的男人了，那么这个课程中应当融入你的生信理念。
1. 这个课程应当是个从零到1的课程，希望你能清空你已经安装的相关生信软件，从找数据，下载数据，安装软件，完全反映出一个生信工作者工作的流程，其中安装软件部分，希望能够融入你准备的软件安装技能，比如conda安装时如何避免Jimmy说的环境污染的问题。
1. 既然你也擅长Linux，那么其中应当有一节讲述linux 的基本技能，或者shell 编程，你也可以选择在处理数据时演示出来。
1. 在课程开始的时候，要讲解ATAC的原理，以及我们要处理的这个数据的来龙去脉（我会阅读文献，帮助处理），要讲解ATAC的应用前景，这个可以参考嘉因生物。
1. 简单说来，希望这是一个实时演示的过程，数据如何下载，每个软件的安装，最好的情况时，一个完全不懂的孩子被讲懂了。
1. 在正文讲完了之后，希望能够演示一下，如何用网页工具或者现成的pipeline来处理数据，网页工具就是你发过的类似的ATAC的帖子，现成的pipeline就是ezATAC那个包。我们希望能够用不同的方法得出相同的结果，如果结果不一样，希望你能给出解释，并且给大家上一课。
1. 在这篇文章处理完了之后，我们就开始上线课程。那么过一段时间之后，我们在后面更新别的数据集的处理方法，一来算是实战复习，而是增加卖点。我们说了的哈，这是你的正式作品，应当认真对待。

## 背景： 染色质和染色体的结构和功能

每一条染色单体由单个线性DNA分子组成。细胞核中的DNA是经过高度有序的包装，否则就是一团乱麻，不利于DNA复制和表达调控。这种有序的状态才能保证基因组的复制和**表达调控**能准确和高效进行。

包装分为多个水平，核小体核心颗粒(nucleosome core particle)、染色小体(chromatosome)、 30 nm水平染色质纤丝(30 nm fibre)和高于30 nm水平的染色体包装。在细胞周期的不同时期，DNA的浓缩程度不同，间期表现为染色质具有转录活性，而中期染色体是转录惰性。细胞主要处于分裂间期，所以DNA大部分时间都是染色质而不是染色体，只不过大家喜欢用染色体泛指染色质和染色体。

![DNA packaging](http://oex750gzt.bkt.clouddn.com/17-12-15/75842695.jpg)

很久之前大家喜欢研究**中期的染色体**，原因是光学显微镜只能看的到这种高度浓缩状态的DNA结构。不过**中期染色体**在转录上是惰性的，没有研究**间期染色体**的意义大。后来技术发展了，大家就开始通过荧光蛋白标记技术以及显微镜技术研究间期染色质的三维结构和动态。比如说，**间期染色体**其实并非随机地弥漫在细胞核中，不同的染色体占据相对独立的空间，染色体在细胞核所占的空间称之为**染色体领地**(chromosome territory, CT)。研究发现，贫基因(gene-porr)的染色体领域一般倾向于靠近核膜，而富含基因(gene-rich)的染色体领地通常位于细胞核内部。这也反应了人类社会的情况，富人处于核心区，穷人在边缘地带。

除了染色体细胞核内的三维结构外，还需要谈谈和转录调控相关的染色质的核小体。用**内切核糖酶**--微球菌核酸酶(micrococcal nuclease, **MNase**, MN酶)处理染色质可以得到单个核小体。**核小体**是染色质的基本结构，由DNA、蛋白质和RNA组成的一种致密结构。组蛋白是由2个H3-H4二聚体，2个H2A-H2B二聚体形成的八聚体，直径约为10 nm， 组蛋白八聚体和DNA结合在一起形成的核心颗粒包含146bp DNA。DNA暴露在核小体表面使得其能被特定的核酸酶接近并切割。

![核小体](http://oex750gzt.bkt.clouddn.com/17-12-15/86776671.jpg)

染色质结构改变会发生在与转录起始相关或与DNA的某种结构特征相关的特定位点。当染色质用**DNA酶I(DNase)**消化时，第一个效果就是在双链体中特定的**超敏位点(hypersenitive site)**引入缺口，这种敏感性可以反应染色质中DNA的可接近性(accessible)，也就是说这些是染色质中DNA由于未组装成通常核小体结构而特别暴露出的结构。

> 许多超敏位点与基因表达有关。每个活性基因在启动子区域都存在一个超敏位点。大部分超敏位点仅存在于相关基因正在被表达的或正在准备表达的细胞染色中；基因表达不活跃时他们则不出现。

## 染色质开放区域和ATAC-seq

背景已经谈到，超敏位点和基因表达有关，并且超敏位点反应了染色质的可接近性。也就可以反推出“可接近”的染色质结构区域可能与基因表达调控相关。因此，高通量测序发展后，就有文章开发基于能切割出单个核小体的MNase, 能识别超敏位点的DNase的文库建立方法。

但是这些文章要求的细胞数比较多，试验步骤也比较繁琐，直到2015年，一篇文章[Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position](https://www.nature.com/articles/nmeth.2688.pdf)使用了超敏Tn5转座酶切割染色质的开放区域，并且加上接头(adapter)进行高通量测序方法，这就是ATAC-seq(assay for transposase-accessible chromatin using sequencing )

![Tn5转座酶切割](http://oex750gzt.bkt.clouddn.com/17-12-15/27788379.jpg)

DNase-seq, MNase-seq和ATAC-Seq三种方法的测序区域以及短读在基因组上的分布情况见下图：

![不同酶切分析的peak差异](http://oex750gzt.bkt.clouddn.com/17-12-16/85242509.jpg)

图片来源于[Reveling in the Revealed](https://www.the-scientist.com/?articles.view/articleNo/44772/title/Reveling-in-the-Revealed/)

文章证明ATAC-seq和已有DNase-seq以及FAIRE-seq拥有相似的信噪比，与不同来源的DNase-seq也有很高的相关性，且组间可重复性高(R=0.98)。之后通过数据分析，得到了如下结论：

- ATAC-seq insert sizes disclose nucleosome positions: ATAC-seq 插入片段大小揭示了核小体的位置
- ATAC-seq reveals patterns of nucleosome-TF spacing: ATAC-seq揭示了核小体-TF 间隔模式
- ATAC-seq footprints **infer** factor occupancy genome wide：ATAC-seq的印记能推测出全基因组转录因子位置
- ATAC-seq enables epigenomic analysis on clinical timescales: ATAC-seq 提供为临床提供了时间维度表观分析策略

也就是说ATAC-Seq能帮助你从全基因组范围内推测可能的转录因子，还能通过比较不同时间的染色质开放区域解答发育问题。

## 数据分析概要

ATAC-seq数据分析大致分为如下几个步骤：

1. FASTQ数据预处理：去接头、过滤低质量reads(mean phred  quality > 20) 、去除高比例N(>10%)reads
1. 序列比对和比对后排序(Bowtie2 > 2.3.4, SAMtools > 1.7)
1. 比对后过滤： 
   - 移除未比对序列、非主要联配、低质量联配和重复(-F 1804 for PE)
   - 保留正确配对短读(-f 2)
   - 移除多比对位置序列(MAPQ < 30 )
1. bam转换成bigwig用于可视化展示
1. 使用MACS2进行peak calling

## 数据分析实战

### 基本信息

以2018年发表在Nature的文章“Chromatin analysis in human early development reveals epigenetic transition during ZGA” 为例，介绍ATAC-Seq数据分析的套路。这篇文章在原来ATAC-seq的基础上开发了一种miniATAC-seq技术， 主要是优化DNA纯化的步骤，使得只需要最低只需要20个细胞就能建库。

在该技术下，作者做了人类胚胎2细胞、8细胞、人胚囊萌发第五天的内细胞团（ICM)和人胚胎干细胞（ES）的ATAC-seq和mRNA-seq。作者检查了几个已知和胚胎发育有关的基因，如_NANOG_, _POU5F1_, 证明了方法的可靠性。

![人类早期胚胎发育ATAC-seq可视化](http://oex750gzt.bkt.clouddn.com/18-5-22/95881386.jpg)

### 数据下载

本篇文章提供的GE编号为**GSE101571**，包括人类和小鼠不同时期的ATAC-seq和mRNA-seq共计55个样本，数据量非常的大，

### 数据预处理

数据预处理步骤分为：质量控制，原始序列比对，比对后去除重复序列和细胞器序列。当然在这之前，先得做一下准备工作，创建工作环境，从SRA下载数据并进行数据解压。

**质量控制**：在数据分析之前先要大致了解手头数据的质量，目前基本就用fastqc了

![multiqc展示](http://oex750gzt.bkt.clouddn.com/17-12-16/20840134.jpg)

FastQC结果大部分都过关，除了在read的各位置碱基含量图上fail。具体原因我还不知道，文章中并没有提到要对原始数据进行预处理。

**序列比对**: 目前在Peak-calling这里分析的流程中，最常用的比对软件就是Bowtie, 分为Bowtie1和Bowtie2，前者适合25~50bp，后者适用于 > 50bp的情况。此处分析的read长度为50 bp，因此选择bowtie1。

> 使用之前建议先看看两个工具的手册，了解一下参数说明。我也是通过Bowtie2的手册才发现Bowtie2和Bowite其实是两个用途不同的工具，而不是说bowtie2是用来替代bowtie1.

有文献报道使用Bowtie比对到拟南芥参考基因组TAIR10，参数为-X2000, 允许长达2 Kb的片段， -m1仅保留唯一联配。

**比对后去(标记)重复**：如果不是PCR-free的建库方法，会有大量重复的read，需要**标记**或**去除**重复。这一步在比对的时候用`samblaster`共同完成。

初步比对后，可以统计下比对到organellar genomesh和nuclear genome的read数量。这个工具可以在shell脚本中用samtools处理，也可以用python的pysam模块.Pysam封装了htslib C-API，提供了SAM/BAM/VCF/BCF/BED/GFF/GTF/FASTA/FASTQ的操作，最新版本已经支持Python3,强烈推荐学习。

![统计比例](http://oex750gzt.bkt.clouddn.com/17-12-18/27220737.jpg)

可能这是作者第一次采用INTACT的方式从干细胞和叶肉细胞中提取细胞核做ATAC-Seq，由于植物细胞本身原因，大量read比对到了细胞器上。

**过滤细胞器reads**：由于细胞器DNA蛋白结合少，所以显然更容易被Tn5 转座酶切割，普通的ATAC-Seq的read就会有大量是细胞器的DNA，这就是为啥需要用INTACT技术。

**过滤低比对质量reads**：根据需要可以过滤掉一些Mapping Quality过低的reads，目前看到的标注有MQ>1, MQ >10.

### Peak Calling

**bigwig定量文件**: 使用deepTools进行标准化和可视化, 一般以RPKM做标准化，默认bin为50。

**Peak Calling**: Peak Calling所用软件和ChIP-Seq是一致的，目前大多为HOMER和MACS2。

```bash
macs2 callpeak -t  # 处理文件
    -f BAMPE  # 数据类型
    --keep-dup all # 保留PICARD和SAMTOOL的重复标记
    -g 1.2e8 #有效基因组大小，人类为2.7e9，拟南芥为120M
    --outdir
    -n #实验名
    -p # p值
    --broad # 根据相联的附近高度富集区域找broad peak
    --broad-cutoff

```

**数据可视化**: Peak Calling结果的可视化作图

- [deepTools](http://deeptools.readthedocs.io/en/latest/)
- [ngsplot](https://github.com/shenlab-sinai/ngsplot)
- [CEAS](http://liulab.dfci.harvard.edu/CEAS/usermanual.html)

### 下游分析

转录因子(transcription factor, TF)：凡是转录起始过程中必须的蛋白质，只要它不是RNA聚合酶的组分，就可将其定义为转录因子。许多转录因子是通过识别DNA上的顺式作用元件而发挥功能，但结合DNA并不是转录因子的唯一作用方式。转录因子还可通过识别另一种因子发挥作用，或识别RNA聚合酶，或只是和其他几种蛋白质一起组成复合体。

在真核生物中，转录因子需要结合到开放的染色质结构，因此核小体八聚体先在启动子上被去除才能让转录因子进行结合。

- ATAC-Seq peak-calling: ZINBA 参数 window size=300bp, offset=75bp. 保留先验概率大于0.8.
- 基于染色质注释的ATAC-Seq插入大小富集分析：计算和各个染色质状态重叠的PE序列片段长度分布。
- 核小体位置：首先将read分组， 低于100bp的read为零核小体； 180 ~ 247 bp单核小体；315 ~ 473为双核小体；558 ~ 615bp为三核小体（原因如下图）。双核小体reads被分成2个reads，三核小体被分成3个reads。reads使用Danpos和Dantools(-p 1, -a 1, -d 20, -clonalcut 0)进行分析。使用零核小体作为背景。

![核小体分布](http://oex750gzt.bkt.clouddn.com/17-12-15/9883328.jpg)

## ATAC-Seq分析工具总结

简单罗列一下ATAC-Seq数据分析会用到的工具以及用途

- 数据预处理：Bowtie1,Bowtie2,BWA-mem, SAMtools,sambamba, picard,samblaster
- Peak Calling: HOMER-findPeaks, MACS2, HOTSPOT
- 差异Peak分析： edgeR, HOMER-getDifferentialPeak,HTSeq + DESeq2
- GO分析: GeneCodis 3.0, AgriGO
- KEGG分析
- 根据距离最近的TSS注释Peak： HOMER-annotatePeaks， ChIPseeker, PeakAnnotator
- 聚类分析：GeneCluster 3.0(k-means)
- 数据可视化工具：IGV, DeepTools, Java Treeview,PAVIS web tool
- Motif分析：HOMER-findMotifsGenome.pl 基于HOMER数据库,MEME-ChIP pipeline(基于repeat-masked fasta), DREME, MEME, CentriMo,RSAT(Regulatory Sequence Analysis Tools),MEME-Suite
- 转录因子footprints: pyDNase
- 直系同源基因检测：CoGe SynFind, CoGe SynMap
- 数据库： 已知motif数据库(Cis-BP,DAP-seq)
- 蛋白互作分析： STRING
- 结合位点的转录因子预测：FIMO



**ATAC-seq** 从分析上来看和**ChIP-seq**没啥区别，都是peak-calling，也就是从比对得到BAM文件中找出reads覆盖区，也就是那个峰。那么问题集中在如何找到peak，首先得问问自己peak的到底是啥?

>“Peak不就是HOMER/MACS2这些peak-finder工具找到的结果吗？”

找Peak就好像找美女，你觉得美女要手如柔荑，肤如凝脂，领如蝤蛴，齿如瓠犀，螓首蛾眉，巧笑倩兮，美目盼兮。但实际情况下，是先给你看一个长相平平的人或者有点缺陷的人，然后再把那个人PS一下，你就觉得是一个美女了。理想情况下， peak应该是一个对称的等腰三角形，并且底角要足够的大。实际情况下是稍微不那么平坦似乎就行了。

首先了解一下ChIP-seq分析时会遇到的三种peak模式。第一种是单个转录因子特异性结合基因组某个区域产生的窄峰, 宽度约为100 bp; 第二种是染色质标记在基因组上的结合情况，如组蛋白修饰H3K27me3，用于修饰的核小体的靶向范围可以不同，因此修饰可以是局部事件，例如相对于启动子中核小体，也可以是大范围事件，延伸覆盖一个大的染色质结构域，所以会得到一个宽峰，长度甚至是100 kbp.  第三种则是混合模式，既有窄峰也有宽峰，数量级大概为10 Kbp.

![ChIP-seq三种peak模式](http://oex750gzt.bkt.clouddn.com/18-5-22/56897678.jpg)

那么ATAC-seq会得到什么样的peak模式图呢？



Tn5酶只能够接近开放染色质(open chromatin)区域，该区域的核小体之间的距离可大可小，因此Tn5酶切割染色体后，会有如下几个情况：

- 两个核小体间DNA，序列长度X(X > 39 bp)

- 刚好一个核小体, 序列长度约为147 bp
- 一个核小体加一段间隔区，序列长度约为147 + X bp
- 连续N个核小体，序列长度约为 (147 X N) bp
- 连续N个核小体加一段间隔区，序列长度约为 (147 X N) + X bp

假设目前已经找到了peak，这是不是意味着我们找到转录因子了？不好意思，这不存在的，因为ATAC-Seq只是找到了全基因组范围的开放区域，而这些开放区域的产生未必是转录因子引起，所以需要一些预测性工作。

在一段展开的DNA上，约119 $\AA$ 的转座同型二聚体(灰绿色)会结合到基因组序列上(蓝色)，酶的核心区(41 $\AA$)作用在9 bp 序列上，同时两翼各有将近 10 bp序列(长度约34$\AA$)，这20 bp序列由于空间阻碍不可能被转座酶结合，所以不会被攻击。由于核心区域在这个过程中出现了两次，因此转座事件最小间隔空间大约为38 bp.

![Tn5酶作用方式](http://oex750gzt.bkt.clouddn.com/18-5-22/85042079.jpg)

## 参考文献

- Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNDNA-binding proteins and nucleosome position
- Sequencing depth and coverage: key considerations in genomic analyses
- Chromatin analysis in human early development reveals epigenetic transition during ZGA
- Rapid, low-input, low-bias construction of shotgun fragment libraries by high-density in vitro transposition
- [ENCODE 分析流程](https://www.encodeproject.org/atac-seq/#standards)
- [Cistrome 分析流程](http://cistrome.org/db/#/about)
- [pyflow-ATACseq](https://github.com/crazyhottommy/pyflow-ATACseq)