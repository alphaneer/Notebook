---
title: 生物信息学数据
author: xuzhougeng
tags: unix, 生物信息
notebook: *NIX基础
---
# 生物信息学数据

数据分析最重要的前提就是数据，毕竟巧妇难为无米之炊。目前生物领域产生海量的数据，比如说NCBI上的1000genomes数据就将近10Pb(1Pb=1024Tb)

## 使用wget/curl下载数据

`wget`和`curl`在大部分\*nix系统都有，大部分互联网的资源都可以使用它们进行下载。`wget`比较常用的参数如下

```bash
-A/--accept: 使用*,?,[,]等构建允许下载的模式
-R/--reject: 使用*,?,[,]等构建不允许下载的模式
-nd/--no-directory: 不需要按照远程文件结构下载数据
-r/--recursive: 递归下载
-np/--no-parent: 不要移动到上级文件夹
--limit-rate: 限定下载速度
--usr=usr: 用户
--ask-password: 认证密码
-O: 定义下载后的文件名
```

以从<http://202.127.18.228/RicePanGenome/>下载水稻泛基因组深度测序组装的66个水稻品系的contig为例

```bash
mkdir rice_contigs && cd rice_contigs
wget -r 1 -np -nd -A *.fa.gz http://202.127.18.228/RicePanGenome/
```

![downloads](http://oex750gzt.bkt.clouddn.com/18-1-20/30535265.jpg)

wget能够通过FTP和HTTP协议下载数据，在网页中递归爬取数据。curl和wget有点不一样，它默认下载数据到标准输出，并且支持SFTP和SCP协议。例如从phytozome里下载数据需要登陆，官方推荐的命令行下载方式就是curl

```bash
# 保存cookie
curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=USER_NAME' --data-urlencode 'password=USER_PASSWORD' -c cookies > /dev/null
# 下载用于保存下载的XML文件
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism=PhytozomeV10' -b cookies > files.xml
```

## 从GEO上下载数据

GEO是用来保存芯片数据、高通量测序数据和其他高通量基因组学数据的平台，提供下载数据、上传数据和数据查询和分析服务，这里主要关注如何利用GEO下载数据。

GEO是数据按照三种方式归档，平台(GPLxxx), 样本(GSMxxx)和系列(GSExxx)。文章中最常出现的形式是GSExxx, 例如GSE101571, 归档了文章中所有实验用到的数据的编号。我们的目的是从这个标号查询可供下载的形式。

GSE101571的查询地址在<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101571>，网页主要是方便人类阅读，爬取则依赖更加规则化的文本

![归档](http://oex750gzt.bkt.clouddn.com/18-5-7/32056722.jpg)

- SOFT: Simple Omnibus in Text Format
- MINiML: MIAME Notation in Markup Language, pronunced minimal, 一类XML扩展用于渲染SOFT

从SOFT的下载链接<ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101571/soft/GSE101571_family.soft.gz>，可以找到规则的命名格式，

找到里面的SRX，然后基于SRX用E-Utils找SRR编号<https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=SRX3884851>

<https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=5337710>

## 使用Aspera服务下载数据

```bash
ascp -v -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l200m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRRXXX/SRRXXXYYY/SRRXXXYYY.sra .
```

```bash
ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRRXXX/SRRXXXYYY/SRRXXXYYY_1.fastq.gz .
```