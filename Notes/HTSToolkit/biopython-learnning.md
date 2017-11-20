---
title: biopython的学习笔记
tags: HTSTookit, Python
notebook: 工具笔记
---

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

# biopython的学习笔记

本文基于\*nix系统的Python2，但是代码也能在Python3运行。

## 序列读取

biopython的`Bio.SeqIO`提供了一个统一接口用于读取/输出不同格式的数据, 另一个和它相近的则是`Bio.AlignIO`用于处理序列联配文件。除非你的输入文件格式太过标新立异，甚至都要用到专门的函数，否则都建议`Bio.SeqIO`搞定。

`Bio.SeqIO`目前支持的格式有: abi, ace, clustam embl, fasta, fastq(sanger, solexa, illumina), genbank/gb, ig, imgt, nexus, pdb-seqres, pdb-atom, phd, phylip, pir, seqxml, sff, stockhold, swiss, tab, qual, uniprot-xml. 覆盖了几乎所有你可能遇到的格式(并没有, sam/bam,vcf/bcf文件需要Pysam用于操作).

以拟南芥参考基因组为例(下载方式为`wget -c -4 -q http://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas  &`)介绍使用方法。

```python
from __future__ import print_function
from Bio import SeqIO

```

### 生物序列对象

### 序列注释对