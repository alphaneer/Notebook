---
title: Massive genomic variation and strong selection in Arabidopsis thaliana lines from Sweden
tags: 群体遗传学,拟南芥,自然变异 
notebook: 文献笔记
---
<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

# 出现在瑞典拟南芥中的大规模基因组变异和严格选择

## 基因组测序

超声破碎，peak为400bp。选择450~800 bp的片段，双端测序，长度为76~100 bp.

## 多态性检测

### 初始read比对

#### read比对和SNP探索

比对: BWA(version 0.5.9) 4% mismatch和1 gap。
去重: samtools rmdup, 仅除去文库准备或测序导致的重复，保留高度重复区域(ribosome repeast centromere repeats)
SNP Indel calling: GATK UnifedGenotyper + IndelGeneotyper. GATK realignment调整已有比对中的variants，然后用samtools call SNP.
过滤： samtool.pl varFilter

#### 质量计算和过滤

BWA计算的MQ(mapping-quality)更倾向于捕捉重复而不是错配。保留两个版本，

- 原始版本
- 仅仅保留Q>30

过滤之后的基因组是原有的87%， SNP和indels是原来的86%和85%.

### 基因参考基因组结构变异识别

- 局部重排寻找短indel:
- PE reads 寻找大SV
- 通过覆盖程度寻找copy number variants(CNV): 计算1kb windows的覆盖率，然后根据周围3Mb windows对其标准化。
