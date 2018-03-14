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

## 使用Aspera服务下载数据

```bash
ascp -v -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l200m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRRXXX/SRRXXXYYY/SRRXXXYYY.sra .
```

```bash
ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRRXXX/SRRXXXYYY/SRRXXXYYY_1.fastq.gz .
```