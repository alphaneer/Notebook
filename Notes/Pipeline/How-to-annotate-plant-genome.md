---
title: 如何对基因组序列进行注释
tags: 基因组
notebook: 分析流程
---
# 如何对基因组序列进行注释

基因组组装完成后，或者是完成了草图，就不可避免遇到一个问题，需要对基因组序列进行注释。注释之前首先得构建基因模型，有三种策略：同源预测(homology-based prediction), 从头注释(_de novo_ prediction)和基于转录组预测(transcriptome-based prediction)，然后才是功能注释，蛋白功能域注释，基因本体论注释，通路注释

## 几个例子

### 陆地棉基因组注释

文章标题为“Sequencing of allotetraploid cotton (Gossypium hirsutum L. acc. TM-1) provides a resource for fiber improvement”.

**同源注释**：从Phytozome上下载了7个植物的基因组蛋白序列(Arabidopsis thaliana, Carica papaya, Glycine max, G. raimondii, Populus trichocarpa, Theobroma cacao and Vitis vinifera), 使用 _TblastN_ 将蛋白序列比对到组装序列上，E-value的阈值为1e-5. 将不同蛋白的BLAST的hits用 _Solar_ 软件进行合并。_GeneWise_ 根据每个BLAST hit的对应基因区域预测完整的基因结构。

**从头预测**：先得构建repeat-mask genome， 在这个基础上就用 _August_, _Genescan_, _GlimmerHMM_, _Geneid_ 和 _SNAP_ 预测编码区

**转录组预测**：用Tophat将RNA-seq数据比对到组装序列上，然后用cufflinks组装转录本形成基因模型。

综上，使用 _EvidenceModeler(EVM)_ 将上面的结果组装成非冗余的基因结构。进一步根据Cscore > 0.5，peptide coverage > 0.5 和CDS overlaping with TE进行筛选。还有过滤掉超过30%编码区被Pfam或Interprot TE domain的注释的基因模型。

这些基因模型使用BLASTP进行功能注释，所用数据库为SWiss-Prot和TrEMBL.蛋白功能使用InterProScan和HMMER注释，数据库为InterPro和Pfam。GO注释则是直接雇佣InterPro和Pfam注释得到的对应entry。通路注释使用KEGG数据库。

### Cardamine hirsuta基因组注释

文章标题为“The Cardamine hirsuta genome offers insight into the evolution of morphological diversity”。

**同源注释**：使用 _Genomethreader_ 以拟南芥为剪切模型，以及PlantsGDB resource上 _Brassica rapa_ (v1.1), _A. thaliana_(TAIR10), _A. lyrata_ (v6), _tomato_ (v3.6), _poplar_ (v2) 和 _A. thaliana_ (version PUT-169), _B. napus_ (version PUT-172) EST assemblies 的完整的代表性蛋白集。

**转录本预测**： 将 _C. hirsuta_ RNA-seq数据比对到基因序列，然后用cufflinks拼接

**从头预测**：转录本预测得到的潜在蛋白编码转录本使用网页工具 _ORFpredictor_ 进行预测， 同时用 _blastx_ 和 _A. thalina_ 进行比较，选择90%序列相似度和最高5%长度差异的部分从而保证保留完整的编码框(有启动子和终止子)。 这些基因模型根据相互之间的相似度和重叠度进行聚类，高度相似(>95)从聚类中剔除，保证非冗余训练集。为了训练gene finder, 它们选随机选取了2000个位点，20%是单个外显子基因。从头预测工具为_August_ , _GlimmerHMM_, _Geneid_ 和 _SNAP_ . 此外还用了Fgenesh+, 以双子叶特异矩阵为参数进行预测。

最后使用JIGSAW算法根据以上结果进行训练，随后再次用JIGSAW对每个基因模型计算统计学权重。

可变剪切模型则是基于苗、叶、花和果实的RNA-seq比对组装结果。

GO注释使用AHRD流程(https://github.com/groupschoof/AHRD/)

### 小结

当检索这些工具的时候，我不经意或者说不可避免就遇到了这个网站 <http://www.plantgdb.org/> , 一个整合植物基因组学工具和资源的网站，但是这个网站似乎2年没有更新了。当然下面这个网站更加重要<http://bioservices.usd.edu/gsap.html>, 他给出了一套完整的注释流程以及每一步的输入和输出情况。

2017年在《Briefings in Bioinformatics》发表的"Plant genome and transcriptome annotations: from misconceptions to simple solution" 则是从五个角度对植物基因组注释做了很完整的总结

- 植物科学的常见本体
- 功能注释的常用数据库和资源
- 已注释的植物基因组意味着什么
- 一个自动化注释流程
- 一个参考流程图，用来说明使用公用数据库注释植物基因组/转录组的常规步骤

![注释流程图](http://oex750gzt.bkt.clouddn.com/18-1-23/36791187.jpg)
