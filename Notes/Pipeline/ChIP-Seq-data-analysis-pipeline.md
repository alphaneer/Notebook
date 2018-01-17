---
title: ChIP-Seq数据分析
author: xuzhougeng
tag: ChIP-Seq, epigenetic
notebook: 分析流程
---
# ChIP-Seq数据分析流程

## 简介

利用染色质免疫沉淀(ChIP)研究转录因子在基因组的位置的方法早就10多年前就有。当ChIP与DNA杂交芯片结合时就成了ChIP-chip，而和高通量测序技术结合时就是ChIP-Seq.

目前ChIP-Seq被广泛用于研究不同物种中的**转录因子**、**组蛋白修饰**、**染色质修饰复合体**和其他染色质相关蛋白。下图是ENCODE对不同测序技术以及应用的总体介绍，从中也能看出ChIP-Seq主要对蛋白结合的启动子和大规模调节因子进行测序。

![](http://oex750gzt.bkt.clouddn.com/18-1-9/49684886.jpg)

目前ENCODE和modENCODE已经在4个物种(果蝇、线虫、老鼠和人类)的100多种细胞上做了上千次的ChIP-Seq实验，研究了超过140多不同的转录因子和组蛋白修饰。并且使用的是多个独立的数据结果和分析流程。综上，他们给出了一套实验设计和数据分析指南，重点是如下几个方面

- 免疫沉淀特异性和质量
- DNA测序深度影响
- 数据集的打分和评价
- 合适的对照组、生物学重复和数据报告

## ChIP工作流程

![](http://oex750gzt.bkt.clouddn.com/18-1-9/99351783.jpg)

## 分析过程

### 前期准备

分别下载人类参考基因组序列及其注释信息，并且构建好对应的索引。尽管hg38(对应GRCh38)可能是未来的趋势，但是目前看来传统的惯性太过强大，hg19(对应GRCh37)依旧是主流。

```bash
# build reference
## fasta
cd ~/reference/genome
mkdir hg19 && cd hg19
curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar xf chromFa.tar.gz
mv chromFa.tar.gz ../
## 将解压缩得到单个染色体序列合并成一个，并删除
cat chr*.fa > hg19.fa
rm -f chr*.fa
## annotation
cd ~/reference/gtf
mkdir hg19 && cd  hg19
curl -O ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.annotation.gtf.gz
```

建立索引或者下载索引。由于人类参考基因组较大，从头建立所以比较耗费时间，因此可以直接选择下载bowtie2已经建立的hg19索引文件.

```bash
cd ~/reference/index/bowtie2
mkdir hg19 && cd hg19
curl -O ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip
unzip hg19.zip
rm hg19.zip
```

在[cistrome](http://cistrome.org/db/)挑选你感兴趣的物种和转录因子，跳转到GSM获取更多实验信息和高通量测序数据下载的编号

![](http://oex750gzt.bkt.clouddn.com/18-1-14/55762783.jpg)

下载数据有三种方案，wget/curl, prefetch, aspera, 此次使用aspera服务

```bash
mkdir -p ~/hs-ChIP-Seq-H3K23me3/analysis/00-raw-data
cd ~/hs-ChIP-Seq-H3K23me3/analysis/00-raw-data
ascp -v -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l200m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR309/SRR3098497/SRR3098497.sra .
```

### 序列预处理

SRA格式的需要先用fastq-dump解压

```bash
# 位于00-raw-data目录下
fastq-dump --gzip --split-3 --defline-qual '+' --defline-seq '@$ac-$si/$ri' SRR3098497.sra
```

解压完之后用`trimmomatic`进行质控, 去除低质量序列和接头序列。所分析的数据来自于**Illumina HiSeq 2000**, 因此接头选择TruSeq3。除非仪器是Illumina的**GA2**系列，需要注意使用TruSeq2，默认就用TruSeq3.

```bash
# 位于analysis目录下
mkdir logs #存放软件运行日志
java -Xmx8G -jar ~/biosoft/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 00-raw-data/SRR3098497.fastq.gz 01-clean-data/SRR3098497_clean.fastq.gz ILLUMINACLIP:/path/to/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:5 TRAILING:10 AVGQUAL:20 MINLEN:80 2> logs/trim.log
# 注 /path/to表示到Trimmomatic-0.36绝对路径
```

序列比对这一步使用bowtie2。 为了提高效率，还可通过管道接`samtools view`和`samtools sort`完成格式转换和排序。

```bash
~/biosoft/bowtie2-2.3.4/bowtie2 -p 10 --sensitive -x ~/reference/index/bowtie2/hg19/hg19 -U 01-clean-data/SRR3098497_clean.fastq.gz 2> logs/align.log | samtools view -b -q 10 | samtools sort > 02-read-alignment/SRR3098497.bam &
```

注：此处过滤了MQ低于10的比对。

可视化

- reads在基因组位置分布统计
- reads相对基因位置分布统计

## peak calling

上面几步都是高通量常规套路，而peak calling则是ChIP-seq/MNase-seq/ATAC-seq所特异的分析方式。peak calling比较常用的软件为[HOMER](http://homer.ucsd.edu/homer/)和[MACS2](https://github.com/taoliu/MACS/). 这里使用MACS进行peak calling， 参数和用法见[如何使用MACS进行peak calling](https://www.jianshu.com/p/6a975f0ea65a)

## peak可视化

目前找到peak后常见的可视化类型有：

- peak长度分布柱状图
- 每个peak的测序情况可视化(IGV, sushi)
- peaks相对基因位置分布统计
- 统计peak在各种基因组区域分布情况
- peak与转录起始位点距离分析

## 下游分析

下游分析的套路主要按照文章[Epigenomic profiling of primary gastric adenocarcinoma reveals super-enhancer heterogeneity](http://dx.doi.org/10.1038/ncomms12983)里的分析方法进行。