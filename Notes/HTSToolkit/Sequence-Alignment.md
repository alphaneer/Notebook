---
title: 序列联配
date: 2017/11/26
tags: alignment
notebook: 生物信息学
categories: HTSToolkit
comments: true
---

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

# 序列联配

序列联配是生物信息学最基础的概念，因为大多数数据分析分析策略都需要使用联配得到的信息。

举个简单的例子，假设你手头上有一些片段'THIS','LI','NE','ISALIGNED', 已知他们来自于一个词，那么原来这个词应该是什么样子。

```shell
--ISALIGNED
THIS-------
-----LI----
--------NE-
```

组成原来单词的过程就是序列联配的所做的事情，就是让碎片化的信息根据相似的部分进行排列，从而推测出原来的数据。这个工作看起来简单，但是涉及到上万条碎片的时候，人脑估计就崩溃了，而且还容易出错，所有就开发相应算法进行实现。

## 联配的表示方法

序列联配结果的展示结果有两种，一种是比较可视化适合人类阅读，一种适合机器处理。默认BLAST得到的结果就是第一种展示方式，如下

```shell
ATGCAAATGACAAATAC
||||   |||.||.|
ATGC---TGATAACT--
```

- '-': 表示gap，可以认为是候选的InDel.
- '|': 表示匹配，一摸一样的那种
- '.': 表示错配，可能有碱基发生了突变。

一般而言，在阅读BLAST结果信息时，都是把上面的序列当作参考序列，下面的序列是自己提供的序列。既然序列比对都是相对而言，那么其实存放形式还可以以更加紧缩的形式存放，这就是CIGAR(compact idiosyncratic gapped alignment report,紧凑的异质间隙排列报告)，处于SAM格式的第六列。

所以第二列的信息，用 `CIGAR` 可以表示为: **4M3D3M1X2M1X1M2D** ，解释一下：

- 4 M(matches)
- 3 D(deletions),
- 3 M(matches),
- 1 X(mismatch),
- 2 M(matches),
- 1 X(mismatch),
- 1 M(match),
- 2 D(deletions)

根据这个CIGAR值和参考序列就能重构出我们自己提供的序列

## 如何计算联配得分

假设你的序列是`ATGAA`，和4个目标序列进行了比较，那么哪个目标序列和自己的序列最相近呢？

```shell
  1      2    3     4
ATGAA ATGAA ATGAA AT-GAA
|.|.| |||.| |||.| || |||
ACGCA ATGCA ATGTA ATCGAA
```

为了挑选出最好的"联配"，就需要按照一定规则对联配结果进行打分，比如说匹配得5分，错配扣4分，出现一个gap后扣10分，并且在已有gap上延伸会多口0.5分。按照这个规则，计算上面联配的得分分别是"7,16,16,15".因此你或许认为第二个和第三个可能才是你的目标序列，毕竟看起来也是这样。但是这只是生信分析得到的一个猜测而已。其实没有最好的联配结果，也没有最好联配相对SCORING选择。改变计分规则，最好的联配结果也会改变。

但是根据一些基本信息，还是能够推理出最优选择。因此在联配的时候，请注意挑选合适的得分矩阵：

* [Selecting the Right Similarity-Scoring Matrix](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3848038/) by William Pearson, the author of the FASTA program.

另外你还可以去[ftp://ftp.ncbi.nlm.nih.gov/blast/matrices](ftp://ftp.ncbi.nlm.nih.gov/blast/matrices)查看目前的得分矩阵类型.

拓展阅读：

- [Wikipedia: Sequence Alignment](https://en.wikipedia.org/wiki/Sequence_alignment)
- [Wikipedia: Gap Penalty](https://en.wikipedia.org/wiki/Gap_penalty)

## 序列联配的算法

在后续的操作需要预先安装一些软件和工具

```shell
conda install -c emboss
mkdir -p ~/bin
curl http://data.biostarhandbook.com/align/global-align.sh > ~/bin/global-align.sh
curl http://data.biostarhandbook.com/align/local-align.sh > ~/bin/local-align.sh
chmod +x ~/bin/*-align.sh
```

下面的例子中使用两条神奇的蛋白序列`THISLINE`和`ISALIGNED`。之所以说是神奇，是因为这完全是认为按照语言习惯对蛋白缩写字母进行排序，当然这个例子最早来自于Marketa Zvelebil和Jeremy Baum的《理解生物信息学》。

> 这两个shell脚本实际上是封装了emboss的needle方法，修改了标准输入和标准输出。

### 全局联配：Global alignments

所谓全局联配，就是尽可能保证两条序列的每个碱基都能配对，因此会有比较多的错配和gap，而且**不会**对**序列两端的gap**进行惩罚。

```shell
# 默认
./global-align.sh THISLINE ISALIGNED
# 自定义gap打开的罚分
./global-align.sh THISLINE ISALIGNED -gapopen 7
```

![](http://oex750gzt.bkt.clouddn.com/17-11-26/49467130.jpg)

### 局部联配：Local alignments

局部联配不是尽可能把两条序列进行配对，而是尽可能的找到那些子区域是能最优联配。并且按照得分矩阵，产生分数在阈值之内的比对结果

```shell
./local-align.sh THISLINE ISALIGNED
a                  7 NE      8
                     ||
b                  7 NE      8
```

这个比对是依赖于EBLOSUM62得分矩阵，默认得分低于11都会被淘汰掉。我们可以下载其他得分矩阵进行比较。

```shell
wget -4 ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/NUC.4.4 -q &
wget -4 ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM30 -q &
wget -4 ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM62 -q &
wget -4 ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM90 -q &
```

比较之后发现，在BLOSUM90下得到的序列长度最长

```shell
./local-align.sh THISLINE ISALIGNED -data BLOSUM90
4 SLI-NE
 :|| ||
3 ALIGNE
```

局部联配适用于寻找两个序列最相似的区域。

> 如果两个序列完全相似，那么局部联配的结果和全局联配的结果将会一致。毕竟每一个地方都是最好的，总体而言也是最好的。如果只是为了实现总体最佳，而不顾及每一个碱基的感受，有一些碱基就得为了实现大局而舍弃自己了。

在全局联配和局部联配之间还有一种混合算法，叫做半全局联配(semi-global alignments), 仅仅知道这么一回事就行。

### 联配的窘境

由于没有上帝视角，我们只能以手头的信息去推测未知的部分。序列联配的时候得到的最优结果仅仅是根据经验设计的算法运行后的结果。给你两条序列“CCAAACCCCCCCTCCCCCGCTTC“和”CCAAACCCCCCCCTCCCCCCGCTTC”，

到底是`./global-align.sh CCAAACCCCCCCTCCCCCGCTTC CCAAACCCCCCCCTCCCCCCGCTTC`得到的结果接近真相

```shell
1 CCAAACCCCCCC--TCCCCCGCTTC     23
  ||||||||||||  .||||||||||
1 CCAAACCCCCCCCTCCCCCCGCTTC     25
```

还是`./global-align.sh CCAAACCCCCCCTCCCCCGCTTC CCAAACCCCCCCCTCCCCCCGCTTC -gapopen 9`得到的结果接近真相呢？

```shell
1 CCAAA-CCCCCCCT-CCCCCGCTTC     23
  ||||| |||||||| ||||||||||
1 CCAAACCCCCCCCTCCCCCCGCTTC     25
```

比对未必能发现真实的插入和缺失，因为算法倾向于用最简单的解释去协调不同的序列，以一句话作为总结吧：

> Alignment reliability depends on the information content of the aligned sequence itself. Alignments that include low complexity regions are typically less reliable. Additional analysis is typically needed to confirm these results.

## 多序列联配

如果要不同物种同一个基因做进化分析，那么就需要先进行多序列联配，然后才能用特定的软件确定不同序列之前的进化关系。常用的软件为: mafft, muscle, clusta-omega, t-coffee等。算法核心思想大多数基于序列递增，但是具体实现有所不同。

以植物转录因子STAT家族的同源序列作为例子，数据来自于北京大学开发的planttfdb, 自行去<http://planttfdb.cbi.pku.edu.cn/family.php?fam=STAT>下载序列。我将序列存放在STAT文件夹下。

一共有214条序列

```shell
# 序列所在文件夹
$ seqkit stat STAT/seq.fas
file     format  type     num_seqs  sum_len  min_len  avg_len  max_len
seq.fas  FASTA   Protein       214  152,000      330    710.3    2,301
```

仅仅挑选出前10条序列

```shell
seqkit seq -n STAT/seq.fas | cut -f 1 -d ' '  | head -n 10 > STAT/ids.txt
seqkit grep --pattern-file STAT/ids.txt STAT/seq.fas > STAT/small.fa
```

使用mafft进行多序列联配

```shell
mafft --clustalout STAT/small.fa > STAT/alignment.maf
```

结果可以用less,head,tail查看，仅选取部分用作展示

![](http://oex750gzt.bkt.clouddn.com/17-11-26/11870199.jpg)