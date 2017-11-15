---
title: 1,135 Genomes Reveal the Global Pattern of Polymorphism in Arabidopsis thaliana
tags: 组装
notebook: 文献笔记 
---
# 实验方法

## Variant calling

GMI-GATK 和 MPI-SHORE两套流程寻找变异位点.

GMI-GATK流程：

- BWA: 序列比对到TAIR10
- samtools: 格式转换，去重
- GATK: indel 局部重比对(local realignment)
    - `UnifiedGenotyper` 寻找indels, 然后用`realignerTargetCreator`生成一组`indelRealigner`所需要的intervals
    - 分别独立地用`UnifiedGenotyper`识别SNP和indels，然后用`CombineVariants`进行合并
- TE-locate: 以1000bp为解析度寻找转座子。

MPI-SHORE流程：

- 比对： BWA sampe(v0.6.2) -n 0.1
- snp calling: SHORE consensus, 得到经验矩阵
- 筛选quality>=25, 替换碱基最小等位频率(minimum allele frequency)大于0.9,转换成VCF

上述两个流程都会得到VCF，仅仅保留一致的交集，且质量大于25. 之后用VCF merge tool合并成 Full Genome VCF. 最后用SnpEff对VCF进行注释。

## 全基因组关联分析

## 群体遗传学分析