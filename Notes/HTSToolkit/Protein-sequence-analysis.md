---
title: 蛋白序列分析
date: 2017/10/12
tags:
  - Sequence Analysis
  - Protein
  - Overview
categories:
  - HTSToolkit
comments: true
---

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

* [蛋白序列分析](#蛋白序列分析)
	* [几个概念](#几个概念)
	* [数据库](#数据库)
	* [分析工具](#分析工具)

<!-- /code_chunk_output -->

# 蛋白序列分析

## 几个概念

- **Homolougous Superfamily**:A homologous superfamily is a group of proteins that share a common evolutionary origin, reflected by similarity in their structure. Since superfamily members often display very low similarity at the sequence level, this type of InterPro entry is usually based on a collection of underlying hidden Markov models, rather than a single signature.
- **Family**: A protein family is a group of proteins that share a common evolutionary origin reflected by their related functions, similarities in sequence, or similar primary, secondary or tertiary structure. A match to an InterPro entry of this type indicates membership of a protein family.
- **Domain**: Domains are distinct functional, structural or sequence units that may exist in a variety of biological contexts. A match to an InterPro entry of this type indicates the presence of a domain.
- **Repeat**t: A match to an InterPro entry of this type identifies a short sequence that is typically repeated within a protein.
- **Site**: A match to an InterPro entry of this type indicates a short sequence that contains one or more conserved residues. The type of sites covered by InterPro are active sites, binding sites, post-translational modification sites and conserved sites.

## 数据库

[InterPro](http://www.ebi.ac.uk/interpro/)：蛋白聚类数据库，以家族，功能域，保守位点进行分类。

[Pfam](http://pfam.xfam.org/): 保守蛋白家族和功能域数据库，是InterPro数据库的其中一部分。

[PDBe](http://www.ebi.ac.uk/pdbe/): 生物大分子和复合体的结构数据库, 数据来源于PDB和EMDB。

[UniProt](http://www.uniprot.org/): 蛋白序列和功能注释的综合性数据库。

[UniProt-GOA](https://www.ebi.ac.uk/GOA): 蛋白序列的高质量基因本体论(GO)注释数据库

## 分析工具

有如下几类分析策略：

- 序列比对，知由来
- 结构聚类，找特征
- 信息注释，寻关联

如果是搜索数据库，从严到松软件选择: BLAST, FASTA > BLAT > PSI-BLAST > HMMER.

> 先以已知序列经BLAST/FATA找到近缘候选序列，就能知道序列的共同的功能域，从Pfam找到相应的HMM Profile 用HMMER搜索远缘基因。

多序列比对软件：Clustal Omega, Kalign, MAFFT, T-Coffee, MUSLE

> 多序列比对一般用于构建进化发育树。

官方提供三种方式使用InterProScan：本地版本，网页服务或网页工具

**网页工具**为<http://www.ebi.ac.uk/interpro/search/sequence-search>， 允许提交长达4000个氨基酸序列的数据。这种方法满足基础使用。

**网页服务**则是使用RESTful, 向服务器提供工作（每次不超过30项），中等规模的数据访问。

**本地服务**需要下载大于5.2G的软件，且要求Linux 64位系统， Perl, Python 2.7.x 以及Java8。**使用方法**为`./interproscan.sh`, 根据提示选择子功能。大规模数据分析。
