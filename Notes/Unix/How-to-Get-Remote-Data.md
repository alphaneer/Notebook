---
title: 生物信息学数据
author: xuzhougeng
tags: unix, 生物信息
notebook: *NIX基础
---
# 生物信息学数据

数据分析最重要的前提就是数据，毕竟巧妇难为无米之炊。目前生物领域产生海量的数据，比如说NCBI上的1000genomes数据就将近10Pb(1Pb=1024Tb)

## 使用Aspera服务下载数据

```bash
ascp -v -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l200m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRRXXX/SRRXXXYYY/SRRXXXYYY.sra .
```

```bash
ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRRXXX/SRRXXXYYY/SRRXXXYYY_1.fastq.gz .
```