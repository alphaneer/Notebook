---
title: 如何构建可重复且稳健的生信项目
author: xu zhougeng
tags: pipeline
notebook: 分析流程
---
# 如何构架可重复且稳健的生信项目

为了保证项目能够重复，需要保证三点：

- 参考基因组版本控制
- 软件版本控制
- 项目结构及版本控制

尽管是三点，但是核心思想就是**版本控制**，同一个项目使用同个物种的不同参考基因组，GO/KEGG数据库的不同版本，同个软件的不同版本等都会引起结果的细微差别。

只要这种细尾差被不会导致结论被完全推翻，那么也就是在误差范围之内。但如果像结构生物学家那样，一些细微差别就是晶体有和无的两个结果，没有任何过渡阶段，那么你就得尽量保证所有因素都要一致。

## 建立统一的参考基因组

以人类为例，目前在[UCSC Browser](http://hgdownload.soe.ucsc.edu/downloads.html#human)上查到人类的参考基因组版本有:GRCh38/hg38,GRCh37/hg19,NCBI36/hg18,NCBI35/hg17等多个版本，不过目前用的最多是前两者。

仅仅有参考基因组也不够，还需要配备对应的注释信息。不同人类参考基因组对应不同的注释信息。由于研究的不断深入，注释信息也在不断更新。在[GENCODE](http://www.gencodegenes.org)中，人类的GRCh38的注释信息截至目前(08.2017)已经更新到27版.

在后续序列比对和变异注释过程中，还需要使用不同的软件根据上述参考基因组和对应注释构建索引，例如BWA, Bowtie2, HISAT2等。

为了对参考基因组进行版本控制，我参考Galaxy项目的文件结构，为研究的拟南芥构建了专门的目录结构存放参考信息，大致如下。

```bash
reference/
├── blastdb
├── genome
│   └── TAIR10
│       ├── Athaliana.fa
│       ├── Athaliana.fa.dict
│       └── Athaliana.fa.fai
├── gtf
│   ├── Phytozome
│   └── TAIR
│       └── TAIR10_GFF3_genes.gtf
├── index
│   ├── BWA
│   │   └── TAIR10
│   │       ├── Athaliana.amb
│   │       ├── Athaliana.ann
│   │       ├── Athaliana.bwt
│   │       ├── Athaliana.pac
│   │       └── Athaliana.sa
│   └── HISAT2
│       └── TAIR10
└── SnpEff
```

目前分为5个部分，介绍如下

- blastdb: 存放不同参考基因组BLAST的核酸和氨基酸索引
- genome: 下级目录为不同参考基因组的版本，然后存放对应的FASTA文件以及字典文件(dict)和索引文件(fai)
- gtf: 下级目录为注释信息的来源，如GENCODE, ENSEMBL,TAIR, Araport等
- index: 下级目录为不同注释软件名
- SnpEff: SnpEff的数据库

具体步骤:

```bash
mkdir -p reference/{genome/TAIR10,index/{BWA/TAIR10,HISAT2/TAIR10},blastdb,gtf/{TAIR,Phytozome},SnpEff}
```

## 统计

## 项目结构以版本控制

为了保证项目的可重复性，就需要为每一个项目创建独立的文件夹，通过不同的文件命名实现脚本的可重复运行。

首先在`.bashrc`中加入如下函数

```bash
pj() { mkdir -p "$1"/{reference/{genome,index},script/{Python,R,Shell},analysis}; }
```

```bash
variant_calling/
├── analysis
│   └── 00-raw-data
│       ├── BC.bg.reads1.fq.gz -> /public/plants/backcross_analysis/Col-0/BC.bg.reads1.fq.gz
│       ├── BC.bg.reads2.fq.gz -> /public/plants/backcross_analysis/Col-0/BC.bg.reads2.fq.gz
│       ├── BC.fg.reads1.fq.gz -> /public/plants/backcross_analysis/BCF2/BC.fg.reads1.fq.gz
│       └── BC.fg.reads2.fq.gz -> /public/plants/backcross_analysis/BCF2/BC.fg.reads2.fq.gz
├── references
│   ├── genomes
│   │   └── TAIR10 -> /public/plants/reference/TAIR10
│   ├── gtf -> /public/plants/reference/gtf
│   └── indexes
└── scripts
    ├── Python
    ├── R
    └── Shell

```