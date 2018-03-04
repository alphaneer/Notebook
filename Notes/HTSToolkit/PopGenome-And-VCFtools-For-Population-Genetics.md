---
title: 群体遗传学分析
tags: 群体遗传学 
notebook: 工具笔记
---
# 群体遗传学分析

， 上游分析流程得到的变异信息总需要进一步的数据挖掘才能找到有趣的故事，但是首先我们得知道可以用什么工具进行什么类型的分析。这篇主要介绍stacks/populationsPopGenome,VCFtools和plink这几个工具可以让我们从数据中挖掘哪些信息。

## stacks/populations

之前只用了stacks的`populations`进行snp过滤和数据导出为常用格式，但其实它还可以计算群体遗传统计值如期望/观测杂合率，π(核苷酸多样性)， Fis（个体,individual相对亚群subpopulation的近交系数），还可以逐对比较所有群体计算Fst，用来判断两个群体间的差异，还可以计算基于单倍体的群体遗传统计值，如单倍型多样性， ΦST, and FST’.

**第一步**，获取VCF文件。其实stacks是可以根据sstacks处理后的数据进行分析，但是VCF可以说是目前许多文件格式的中转站，用的也比较广泛，为了方便后续分析和协作分析，最好还是先得到一个VCF吧。

```bash
min_samples=0.80
min_maf=0.05
max_obs_het=0.80
populations -P 03-stacks-analysis/ref-based/ -r $min_samples --min_maf $min_maf \
--max_obs_het $max_obs_het --vcf &> populations_ref_based.oe
populations -P 03-stacks-analysis/de-novo/ -r $min_samples --min_maf $min_maf \
--max_obs_het $max_obs_het --vcf &> populations_de_novo.oe
```

最后会在`03-stacks-analysis/ref-based/`和`03-stacks-analysis/de-novo/`下都有一个`batch_1.vcf`，这就是后续用来分析的基础。分别改名成`de_novo.vcf`和`ref_based.vcf`便于区分。

> populations支持导出的格式非常的多，如fasta, fasta\_strict, vcf, vcf\_haplotypr, structure,phase,fastphase,beagle,beagle\_phased,plink,hzar,phylip,phylip\_var,phylip\_var\_all,treemix

**第二步**，声明用于配对比较的群体. 这一步生成`-M,--popmap`输入文件。比如说比较淡水鱼(pcr,wj)和海水鱼（cs,sj)，

```bash
sed -r 's/\t(cs|sj)/\tocceanic/; s/\t(pcr|wc)/\tfreshwater/;' info/popmap.tsv \
> info/popmap.oceanic_freshwater.tsv
```

**第三步**，计算群体遗传统计值

```bash
popmap=info/popmap.oceanic_freshwater.tsv
populations -V ref_based.vcf -M $popmap -O ./ -p 2 \
    --fstats -k --sigma 100000
```

对于有参考基因，可以使用`-k`计算kernel-smoothed π, FIS, FST, FST', and ΦST，至于有啥意义，我目前还不懂。此外还可计算**重抽样**的几个统计值，对应`--bootstrap_xx`。启用后计算时间会明显增加。最后会得到如下文件

- ref\_based.fst\_occeanic-freshwater.tsv
- ref\_based.fst\_summary.tsv
- ref\_based.haplotypes.tsv
- ref\_based.hapstats.tsv
- ref\_based.phistats.tsv
- ref\_based.phistats\_occeanic-freshwater.tsv
- ref\_based.sumstats.tsv
- ref\_based.sumstats\_summary.tsv

**第四步**，R可视化展示. 结果必须可视化才能便于我们找到数据里有意思的地方。这一步要用到上一步得到的`ref\_based.phistats\_occeanic-freshwater.tsv`

```R
x = read.delim('ref_based.fst_occeanic-freshwater.tsv')
x.vii = subset(x, Chr=='groupVII')
plot(x.vii$BP, x.vii$AMOVA.Fst, 
     pch=3, cex = 0.5,
     xaxt="n",
     xaxs= 'i',
     yaxs = "i",
     xlab="groupVII Position (Mb)",ylab="Fst")
lines(x.vii$BP, x.vii$Smoothed.AMOVA.Fst,col='blue')
axis(1,at=seq(0,ceiling(max(x.vii$BP)/1000000)*1000000,1000000),
     labels = seq(0,ceiling(max(x.vii$BP)/1000000),1))

```

## PopGenome

PopGenome号称是群体遗传学分析的瑞士军刀，能够处理1000/1001基因组项目的全基因组SNP数据和序列联配后的FASTA数据。在统计检验上，它调用Hudson’s MS 和Ewing’s MSMS程序通过联合模拟进行显著性测试。

由于PopGenome是R包，因此`install.packages("PopGenome")`即可以安装。