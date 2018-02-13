---
title: 对BAM文件按照基因组特征计数
tags: 工具笔记, RNA-seq
notebook: 工具笔记
---
# 对BAM文件按照基因组特征计数

在测序数据比对到参考基因组后,还需要进一步将reads归到对应的基因组特征(genomic features)上。这个步骤通常被称为"read summarization"或是"read quantification"(暂且翻译为短读定量吧). 对于下游游分析不可或缺的一步，无论是基因表达量分析，还是组蛋白修饰分析。这一步会得到每个计数表，记录着每个基因组特征相应的数量。这一步主要的难题是，对于二代测序而言read的长度大多是100bp左右，而一个基因的长度大多是1kb左右，因此当一个read在可变剪接的基因的外显子发生重叠时，read应该如何归类。

这一步定量可以分为三个水平：基因水平(gene-level), 转录本水平(transcript-level)和外显子使用水平(exon-usage-level)。

## 基因水平定量

**基因水平**的计数是最简单的途径，常用工具是**HTSeq-count**和**Subread/featureCounts**，而后者在速度上远远快于前者，结果上不存在太差区别，因此更加推荐使用。**featureCounts**有如下特征：

- 由于考虑到短读中插入缺失,接合点(junction)和结构变异, 能够准确和精确地对短读进行分配
- 运行速度快，官方宣称半分钟处理200万条短读(线程数-T)
- 支持GTF和SAF格式的注释(-a 注释文件)
- 支持链特异性短读定量(-s 0/1/2, 0表示为非链特异性，1表示正链，2表示负链)
- 可以在特征(例如exon)和元特征(例如gene)层面上对基因进行定量(-t exon)
- 在处理存在多个匹配和重叠的短读上高度灵活，可以选择剔除(官方建议，也是默认参数)，可以选择全部保留(-M | -O)或部分保留(-M --fraction | -O --fraction)。全部保留是将每个位置都视为**1**，部分则是**1/x|1/y**, x就是匹配数，y是重叠特征数。
- 可以完全掌控双端测序短读的定量，例如检查是否两端都在同一个基因内
- 通过搜索双端测序的短读在重叠特征(overlap feature)上位置，降低了定量的模糊性(-p -B)
- 允许使用者决定是否对嵌合片段(chimeric fragments)定量(-p -C)
- 自动检测输入格式是BAM还是SAM
- 自动对双端测序数据排序

`featureCounts`的参数非常多，上面已经提及其中几种。对于常规的分析，也不需要了解这些参数，使用默认参数即可。下面是常用的搭配方法：

```bash
# 单端, SAM, 5个进程
featureCounts -T 5 -a annotation.gtf -t exon -g gene_id -o counts.txt mapping_results_SE.sam
# 单端, BAM
featureCounts -a annotation.gtf -t exon -g gene_id -o counts.txt mapping_results_SE.bam
# 多个BAM文件同时处理
featureCounts -a annotation.gtf -t exon -g gene_id -o counts.txt mapping_results1.bam mapping_results2.bam
# 双端，按照fragments而不是reads定量
featureCounts -p -a annotation.gtf -t exon -g gene_id -o counts.txt mapping_results_PE.bam
# 双端，对fragments过滤
featureCounts -p -P -d 50 -D 600 -a annotation.gtf -t exon -g gene_id -o coutns.txt mapping_results_PE.bam
# 不考虑fragments长度，只要fragments都能比对即可
featureCounts -p -B -a annotation.gtf -t exon -g gene_id -o counts.txt mapping_results_PE.bam
# 剔除嵌合片段
featureCounts -p -C -a annotation.gtf -t exon -g gene_id -o counts.txt mapping_results_PE.bam
```

## 转录本水平定量

转录本水平上定量分为两类工具，一类是基于基因组或转录组比对结果，如`stringtie`或`eXpress`，另一类则是不需要比对直接根据转录组定量，如`Sailfish`,`Salmon`和`kallisto`.这些软件要处理的难题就时转录本亚型（isoforms）之间通常是有重叠的，当二代测序读长低于转录本长度时，如何进行区分？这些工具核心算法大多是极大似然法，也就是根据统计学的方法判断这条read最有可能归谁。