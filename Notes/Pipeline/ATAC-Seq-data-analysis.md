---
title: ATAC-Seq数据分析
author: xuzhougeng
tag: ATAC-Seq, epigenetic
notebook: 分析流程
---
# ATAC-Seq 数据分析

## 背景： 染色质和染色体的结构和功能

每一条染色单体由单个线性DNA分子组成。细胞核中的DNA是经过高度有序的包装，否则就是一团乱麻，不利于DNA复制和表达调控。这种有序的状态才能保证基因组的复制和**表达调控**能准确和高效进行。

包装分为多个水平，核小体核心颗粒(nucleosome core particle)、染色小体(chromatosome)、 30 nm水平染色质纤丝(30 nm fibre)和高于30 nm水平的染色体包装。在细胞周期的不同时期，DNA的浓缩程度不同，间期表现为染色质具有转录活性，而中期染色体是转录惰性。细胞主要处于分裂间期，所以DNA大部分时间都是染色质而不是染色体，只不过大家喜欢用染色体泛指染色质和染色体。

![DNA packaging](http://oex750gzt.bkt.clouddn.com/17-12-15/75842695.jpg)

很久之前大家喜欢研究**中期的染色体**，原因是光学显微镜只能看的到这种高度浓缩状态的DNA结构。不过**中期染色体**在转录上是惰性的，没有研究**间期染色体**的意义大。后来技术发展了，大家就开始通过荧光蛋白标记技术以及显微镜技术研究间期染色质的三维结构和动态。比如说，**间期染色体**其实并非随机地弥漫在细胞核中，不同的染色体占据相对独立的空间，染色体在细胞核所占的空间称之为**染色体领地**(chromosome territory, CT)。研究发现，贫基因(gene-porr)的染色体领域一般倾向于靠近核膜，而富含基因(gene-rich)的染色体领地通常位于细胞核内部。这也反应了人类社会的情况，富人处于核心区，穷人在边缘地带。

除了染色体细胞核内的三维结构外，还需要谈谈和转录调控相关的染色质的核小体。用**内切核糖酶**--微球菌核酸酶(micrococcal nuclease, **MNase**, MN酶)处理染色质可以得到单个核小体。**核小体**是染色质的基本结构，由DNA、蛋白质和RNA组成的一种致密结构。组蛋白是由2个H3-H4二聚体，2个H2A-H2B二聚体形成的八聚体，直径约为10 nm， 组蛋白八聚体和DNA结合在一起形成的核心颗粒包含146bp DNA。DNA暴露在核小体表面使得其能被特定的核酸酶接近并切割。

![](http://oex750gzt.bkt.clouddn.com/17-12-15/86776671.jpg)

染色质结构改变会发生在与转录起始相关或与DNA的某种结构特征相关的特定位点。当染色质用**DNA酶I(DNase)**消化时，第一个效果就是在双链体中特定的**超敏位点(hypersenitive site)**引入缺口，这种敏感性可以反应染色质中DNA的可及性(accessible)，也就是说这些是染色质中DNA由于未组装成通常核小体结构而特别暴露出的结构。

> 许多超敏位点与基因表达有关。每个活性基因在启动子区域都存在一个超敏位点。大部分超敏位点仅存在于相关基因正在被表达的或正在准备表达的细胞染色中；基因表达不活跃时他们则不出现。

## 染色质开放区域和ATAC-Seq

背景已经谈到，超敏位点和基因表达有关，并且超敏位点反应了染色质的可及性。也就可以反推出“可及性”的染色质结构区域可能与基因表达调控相关。于是2015年的一篇文章[Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNDNA-binding proteins and nucleosome position](https://www.nature.com/articles/nmeth.2688.pdf)就使用了超敏Tn5转座酶切割染色质的开放区域，并且加上接头(adapter)进行高通量测序。

![](http://oex750gzt.bkt.clouddn.com/17-12-15/27788379.jpg)

![](http://oex750gzt.bkt.clouddn.com/17-12-16/10917527.jpg)

那篇文章通过ATAC-Seq得到了如下结论：

- ATAC-seq insert sizes disclose nucleosome positions
- ATAC-seq reveals patterns of nucleosome-TF spacing
- ATAC-seq footprints **infer** factor occupancy genome wide
- ATAC-seq enables epigenomic analysis on clinical timescales

也就是说ATAC-Seq能帮助你从全基因组范围内推测可能的转录因子，还能通过比较不同时间的染色质开放区域解答发育问题。

## 数据分析概要

在前面的铺垫工作中，一共提到了三种酶,能切割出单个核小体的MNase, 能识别超敏位点的DNase 和ATAC-Seq所需要的Tn5 transposase，这三种酶的异同如下图：

![不同酶切分析的peak差异](http://oex750gzt.bkt.clouddn.com/17-12-16/85242509.jpg)

图片来源于[Reveling in the Revealed](https://www.the-scientist.com/?articles.view/articleNo/44772/title/Reveling-in-the-Revealed/)

分析ATAC-Seq从本质上来看和分析ChIP-Seq没啥区别，都是peak-calling，也就是从比对得到BAM文件中找出reads覆盖区，也就是那个峰。(尴尬的是，这句话对于老司机而言是废话，对于新手而言则是他们连ChIP-Seq都不知道)那么问题集中在如何找到peak，peak的定义是啥?

>“Peak不就是HOMER/MACS2/ZINBA这些peak-finder工具找到的结果吗？”

找Peak就好像找美女，你觉得美女要手如柔荑，肤如凝脂，领如蝤蛴，齿如瓠犀，螓首蛾眉，巧笑倩兮，美目盼兮。但实际情况下，是先给你看一个长相平平的人或者有点缺陷的人，然后再把那个人PS一下，你就觉得是一个美女了。理想情况下， peak应该是一个对称的等腰三角形，并且底角要足够的大。实际情况下是稍微不那么平坦似乎就行了。

假设目前已经找到了peak，这是不是意味着我们找到转录因子了？不好意思，这不存在的，因为ATAC-Seq只是找到了全基因组范围的开放区域，而这些开放区域的产生未必是转录因子引起，所以需要一些预测性工作。

## 数据分析实战

### 基本信息

以目前预发表在bioRxiv的文章“Chromatin accessibility changes between Arabidopsis stem cells and mesophyll cells illuminate cell type-specific transcription factor networks” 为例，介绍ATAC-Seq数据分析的套路。

**GEO编号**：GSE101940，一共6个样本，SRR为SRR5874657~SRR587462

**实验设计**：用INTACT方法提取植物干细胞(stem cell)和叶肉细胞(mesophyll cells)的细胞核，然后通过ATAC-Seq比较两者在转录因子上的差别

**分析流程**：分为数据预处理,Peak-calling和后续分析三步。

- 序列比对：Bowtie2
- 比对后处理： SAMtools和Picard
- BAM可视化： IGV, deeptools
- Peak-Calling: HOMMER
- THS(转座酶超敏位点) 基因区分布： PVIS网页工具，上游2kb 下游1kb
- THS 差异分析： HTSeq + DESeq2
- 转录因子motif分析：RSAT(Regulatory Sequence Analysis Tools), MEME-ChIP, DREME, MEME, CentriMo. 额外数据Cis-BP, DAP-seq
- 调控基因预测：PeakAnnotator
- 蛋白互作分析： STRING
- 结合位点的转录因子预测：FIMO
- 基因富集分析

### 数据预处理

数据预处理步骤分为：质量控制，原始序列比对，比对后去除重复序列和细胞器序列。当然在这之前，先得做一下准备工作，创建工作环境，从SRA下载数据并进行数据解压。

```bash
# 创建项目文件下
mkdir -p ATAC-Seq/{data/raw_data,analysis,script,ref}
# 使用sra-tool prefetch下载数据, 数据保存在~/ncbi/public/sra
for i in `seq 57 62`;
do
    prefetch SRR5874${i} &
done
# 数据下载完，用fastq-dump解压
for i in `seq 57 62`;
do
fastq-dump --split-3 --defline-qual '+' --defline-seq '@$ac-$si/$ri length=$rl' --gzip SRR5874${i} -O data/raw_data &
done
```

**质量控制**：在数据分析之前先要大致了解手头数据的质量，目前基本就用fastqc了

```bash
mkdir -p analysis/fastqc
for i in `seq 57 62`; do fastqc data/raw_data/SRR58746${i}_{1,2}.fastq.gz -o analysis/qc & done
# multiqc汇总
multiqc analysis/qc/ -o analysis/qc/
```

FastQC结果大部分都过关，除了在read的各位置碱基含量图上fail。具体原因我还不知道，文章中并没有提到要对原始数据进行预处理。

![](http://oex750gzt.bkt.clouddn.com/17-12-16/20840134.jpg)

**序列比对**: 目前在Peak-calling这里分析的流程中，最常用的比对软件就是Bowtie, 分为Bowtie1和Bowtie2，前者适合25~50bp，后者适用于 > 50bp的情况。此处分析的read长度为50 bp，因此选择bowtie1。

> 使用之前建议先看看两个工具的手册，了解一下参数说明。我也是通过Bowtie2的手册才发现Bowtie2和Bowite其实是两个用途不同的工具，而不是说bowtie2是用来替代bowtie1.

```bash
# 下载或者链接已有的参考基因组到ref文件下
ln -s ~/db/Genomes/Athalina/TAIR10/Sequence/TAIR10.fa ref/
# 建立BOWTIE2 index，或者下载已有的index, 或软连接已有的索引
bowtie-build --threads 8 ref/TAIR10.fa  ref/TAIR10
# 注意不要用bowtie2-build
# 序列比对
for i in `seq 57 62`;
do
bowtie -p 10 -S -m 1 -X 2000  --sam-RG "ID:sample_${i}" \
--sam-RG "PL:illumina" --sam-RG "SM:SRR58746${i}" \
ref/TAIR10 \
-1 <(zcat data/raw_data/SRR5874658_1.fastq.gz) \
-2 <(zcat data/raw_data/SRR5874658_2.fastq.gz) | \
samtools sort -@ 6 -m 1G -o analysis/BAM/SRR58746${i}_sorted.bam ;
done &
```

这里使用Bowtie比对到拟南芥参考基因组TAIR10，参数为-X2000, 允许长达2 Kb的片段， -m1仅保留唯一联配。

**比对后去重复和细胞器reads**：

### Peak Calling

### 下游分析

- ATAC-Seq peak-calling: ZINBA 参数 window size=300bp, offset=75bp. 保留先验概率大于0.8.
- 基于染色质注释的ATAC-Seq插入大小富集分析：计算和各个染色质状态重叠的PE序列片段长度分布。
- 核小体位置：首先将read分组， 低于100bp的read为零核小体； 180 ~ 247 bp单核小体；315~473为双核小体；558~615bp为三核小体（原因如下图）。双核小体reads被分成2个reads，三核小体被分成3个reads。reads使用Danpos和Dantools(-p 1, -a 1, -d 20, -clonalcut 0)进行分析。使用零核小体作为背景。

![](http://oex750gzt.bkt.clouddn.com/17-12-15/9883328.jpg)
