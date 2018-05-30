# 使用LUMPY检测结构变异

LUMPY是一款基于概率框架检测结构变异(structure variants)的软件, 它根据read-pair, split-read, read-depth和其他先验知识寻找基因组上可能的结构变异。

软件在编译的时候会先安装HTSLIB，根据Makefile, 需要预先安装好curl和zlib. 此外还推荐安装Python的Pysam和Numpy，Samtools(0.1.18+)，SAMBLASTER(0.1.19+), sambamba。安装完成之后可以用测试数据进行测试, 数据下载地址为<http://layerlab.org/lumpy/data.tar.gz>, 要存放在`/lumpy-sv/data`下

`lumpy`基于paired-end reads比对后得到的三类信息推断SV，局部异常的测序深度，不一致(discordant)的联配和断裂的联配(split-read alignment)。局部异常的测序深度比较容易理解，平均30X测序的地方，如果深度大于100X，意味着存在着拷贝数变异，如果深度程度非常低，可能意味着这里存在 **大片段缺失**。不一致的联配和断裂的联配能够提供的信息更多，如果基因组一个区域齐刷刷的截断(如下图)，就意味着这个区域可能存在插入/缺失。当然也有其他可能，当两个read在不同链或者不同染色体时，可能是易位或倒置。

![read截断](http://oex750gzt.bkt.clouddn.com/18-5-25/16719762.jpg)

因此，运行`lumpy`需要预先整理出不一致的短读以及断裂的联配, 如下是`lumpy`提供的数据预处理方式

```bash
# Align the data
bwa mem -R "@RG\tID:id\tSM:sample\tLB:lib" human_g1k_v37.fasta sample.1.fq sample.2.fq \
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view -S -b - \
    > sample.bam

# Extract the discordant paired-end alignments.
samtools view -b -F 1294 sample.bam > sample.discordants.unsorted.bam

# Extract the split-read alignments
samtools view -h sample.bam \
    | scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > sample.splitters.unsorted.bam

# Sort both alignments
samtools sort sample.discordants.unsorted.bam sample.discordants
samtools sort sample.splitters.unsorted.bam sample.splitters
```

让人感兴趣的是`-F 1294`用来提取不一致的联配，用`samtools flags 1294`可以发现1294表示"PROPER_PAIR,UNMAP,MUNMAP,SECONDARY,DUP"，带上`-F`意味着以上这些标记在我们筛选的联配记录中都不会出现，也就意味着筛选的记录要符合下面要求

- 不能是PROPER\_PAIR: 就是比对工具认为都正确比对到基因组上，在同一条染色体，在同一条链的情况，常见的就是83,147和99,163
- 不能同时是UNMAP和MUNMAP，也就是配对的短读至少有一个能够比对到参考基因组上
- 也不能是SECONDARY， 也就是他必须是主要联配
- 光学重复，DUP, 就更加不能要了

于是，经过上一步，那就得到了包含所有数据的sample.bam，不一致的联配sample.discordants.bam 和断裂联配sample.splitters.bam, 使用作者封装好的调用函数进行结构变异检测。

```bash
lumpyexpress \
    -B sample.bam \
    -S sample.splitters.bam \
    -D sample.discordants.bam \
    -o sample.vcf
```

在得到的结构变异基础上，作者推荐是用[SVTyper](https://github.com/cc2qe/svtyper)进行基因型确定。

提高准确度的方法：根据先验剔除已知低复杂区域和高覆盖的区域，见参考资料的lumpy-sv教程。

## 参考资料

- [samblaster关于不一致短读筛选的解答](https://github.com/GregoryFaust/samblaster/issues/33)
- [lumpy-sv教程](https://github.com/arq5x/lumpy-sv)