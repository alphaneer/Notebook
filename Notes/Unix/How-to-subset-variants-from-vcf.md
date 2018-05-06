---
title: 提取vcf文件中一个contig
author: xuzhougeng
tags: unix, 生物信息
notebook: *NIX基础
---
# 提取vcf文件中一个contig

这是一个很小众的需求。大部分变异检测都是基于组装质量比较高的基因组，而不是那种初步拼接的contig。

由于初步拼接的参考序列通常会有成千上万个contig序列，也就导致在VCF的头文件的`##contig=<ID=xxx,length=xxx>`部分会有成千上万个contig，将这个文件加载到IGV时, IGV会去解析VCF，这将会是非常缓慢的过程，最好的策略就是只提取其中一个contig查看。

一开始，我想到的是这个方法，直接用bcftools去提取目标contig

```bash
bcftools view some_variants.bcf contig_1
```

但是，你会发现HEADER部分还是会保留原来的信息。怎么办呢？

我们这里会用到一个Shell高级技巧--subshell。举个例子，比如说同时查看一个文件的头10行和后10行.

```bash
(head -n 10 ; tail -n10 ) < ~/referece/annotaiton/TAIR10/TAIR10_GFF3_genes.bed
```

这里的`()`就是一个subshell, 它的作用就是同时执行`()`中的所有命令。

用到此处，就是同时完成提取VCF的header部分并只保留目标contig和提取目标区间的所有记录这两个任务，然后保存为vcf.gz格式。

```bash
(bcftools view -h some_variants.bcf | grep -v "##contig" | sed "/^##reference/a $(bcftools view -h some_variants.bcf | grep -w "ID=contig_1")" ; bcftools view -H some_variants.bcf contig_1 ) | bcftools view -Oz -o contig_1.vcf.gz
```

上面是一个非常复杂的shell命令，希望你能看懂。

这里面重复出现了两个变量，一个是`some_variants.bcf`，一个是`contig_1`, 我们可以将其用`$1`和`$2`代替，然后添加到你的`~/.bashrc`中，如下所示

```bash
# function
subvcf()
{
    (bcftools view -h $1 | grep -v "##contig" | sed "/^##reference/a $(bcftools view -h $1| grep -w "ID=$2")" ;  bcftools view -H $1 $2 ) | bcftools view -Oz -o $2.vcf.gz && bcftoos index $2.vcf.gz
}
```

最后用`source ~/.bashrc`添加刚才的函数到环境变量中，最后就是`subvcf some_variants.bcf contig_1` 代替之前的长串命令。