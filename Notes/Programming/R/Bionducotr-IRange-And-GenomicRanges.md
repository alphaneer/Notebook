---
title: Bioconductor的地基--IRanges
tags: R, bioconductor
notebook: 工具笔记
---
# Bioconductor的地基--IRanges

Bioconductor是一个开源项目，包括许多R生物信息学包。这里，首先介绍Bioconductor的核心包：

* GenomicRanges: 用于表示和使用基因组范围，genomic ranges
* GenomicFeatures: 用于表示和使用基因模型和基因组其他特性(genes, exons,UTRs,transcripts等）的范围
* Biostrings和BSgenomes:用于在R中操作基因组序列（比如有些包可以从ranges中提取序列）
* rtrcklayer: 用于读取常见的生物信息数据格式，BED，GTF/GFF，和WIG

## 安装和帮助说明

source("http://biocondutor.org/biocLite.R")
biocLite("GenomicRanges")
验证或升级版本：biocValid()、biocUpdatePackages()

## 使用IRanges储存通用数据

IRanges构建语法：`IRanges(start=NULL, end=NULL, width=NULL, names=NULL)`

```r
library(IRanges)
rng <- IRanges(start=4, end=13)
x <- IRanges(start=c(4,7,2,20), end=c(13,7,5,23))
names(x) <- letters[1:4] #letter内置常量，小写字母，LETTER，大写字母，
class(x) # 查看类
#一些常用方法: start(), end(), width(),range()
# 算术、切片、逻辑操作
start(x) + 4
x[2:4]
x['a']
x[start(x) <4]
# merge
a <- IRanges(start=7, width=4)
b <- IRanges(start=2, end=5)
c(a,b)
```

## 基本Range操作，算术、变形、集合

IRanges对象可以通过算术运算增大或缩小范围

```r
x <- IRanges(start=c(40,80), end=c(67,114)
x + 4L
x - 10 L
# 这些运算从两端同时进行
# 截取部分区域，restrict()
y <- IRanges(start=c(4,6,10,12), width=13)
restrict(y,5,10)
# flank()可以提取每个range的两端部分，
flank(x,width=7)
# reduce() 相当于read组装
set.seed(0)
alns <- IRanges(start=sample(seq_len(50),20),width=5) #sample随机取样，seq_len产生范围数据
head(alns,4)
reduce(alns)
# gaps找出不同序列间隔，默认不关注the beginning of the sequences to the start position of the first ranges, same as end
gaps(alns)
# 集合操作:交集intersect, 差级setdiff, 合集union， 逐对操作，pintersect,psetdiff,punion,pgap
```

## 寻找overlapping ranges

寻找overlap是许多基因组学分析任务的必要环节。RNA-seq，overlaps用于分析细胞活动
我们使用findOverlaps()函数寻找两个IRanges对象的overlaps.

```r
qry <- IRanges(start=c(1,26,19,11,21,7), end=c(16,30,19,15,24,8),names=letters[1:6])
sbj <- IRanges(start=c(1,19,10),end=c(5,29,16),names=letters[24:26])
hts <- findOverlaps(qry, sbj)
hts
Show in New WindowClear OutputExpand/Collapse Output
Hits object with 6 hits and 0 metadata columns:
      queryHits subjectHits
      <integer>   <integer>
  [1]         1           1
  [2]         1           3
  [3]         2           2
  [4]         3           2
  [5]         4           3
  [6]         5           2
  -------
  queryLength: 6 / subjectLength: 3
```

Overlaps表示的是query和subject之间的映射关系。每个query根据我们寻找overlaps的方式不同，在不同的subjects可以有多个hits.单个subjects也可以有多个query的hits.举findOverlaps给出的结果[1] 1 1为例，表示为query的1和subject的1存在映射。我们可以使用accessor函数--queryHits(),和subjectHits()查看索引。

```r
names(qry)[queryHits(hts)]  # 分解成index <- queryHts(hts); name <- names(qry); name[index]
names(suj)[subjectHits(hts)]
```

默认是当range任意部分与subject range有重叠，就被认为是overlap，这个默认type为"any"。也可以修改为type=within

```r
hts_within <- findOverlaps(qry, sbj, type="within")
```

within和any是最常用的两种类型，其他的可以通过help(findOverlaps)了解。还有一个参数为select,用来处理多对多的映射关系，默认是select="all",first表示第一个，last表示最后一个，arbitrary表示是任意。
**时刻留意Overlaps，因为有可能会造成结果的巨大差异**
计算overlaps需要强大的计算能力，需要对query和subject进行一一对比。比较明智的方法是使用*interval tree*进行恰当的排序，使用help(IntervalTree)查看帮助，似乎这个功能好像取消了。
运行findOverlaps，我们需要提取Hits对象中的信息,例如之前用到的subjectHits。

```r
as.matrix(hts) # 直接强迫转换成matrix
countQueryHits(hts) # 计算每个query的命中数
setNames(countQueryHits(hts),names(qry)) # 给对象命名
countSubjectHits(hts)
setNames(countSubjectHits(hts),names(sbj)
ranges(hts,qry,sbj)  # 找出qry,sbj的共同之处
```

ranges所创建的对象可以使用R的所有向量数据分析方法。例如summary.
subsetByOverlaps()仅保留的query的子集，countOverlaps()则计算重叠数。

## 寻找邻近Ranges以及计算距离

有三种方法可以完成邻近ranges的寻找:neatest(), precede(), fellow()

```r
qry <- IRanges(start = 6, end = 13, names = 'query')
sbj <- IRanges(start = c(2,4,18,19), end=c(4,5,21,24),names=1:4)
nearest(qry, sbj) #最近的距离
precede(qry,sbj) #寻找前面的ranges
follow(qry,sbj) #寻找后面的ranges
```

这些操作允许vertorization.·qry2 <- IRanges(start=c(6,7), width=3); nearest(qry2,sbj)`.
使用distanceToNearest()和distance()确定两个ranges之间的距离

```r
qry <- IRanges(sample(seq_len(1000),5),width=10)
sbj <- IRanges(sample(seq_len(1000),5),width=10)
distanceToNearest(qry,sbj)
distance(qry,sbj) #仅返回距离
```

## Run Lenth Enconding and Views

IRanges可以操作任意类型的序列的ranges。如果是基因组数据，那么ranges的坐标依赖于核酸和特定的染色体。当然爱有许多其他类型的基因组数据，例如:

* Coverage, 一定长度序列中ranges的重叠深度
* Conservation tracks, 两个物种间，碱基与碱基之间的进化保守计分
* Per-base pair estimates of population genomics summary statistics like nucleotide diversity.

### Run-length encoding and coverage()

序列越长所占的容量越大，为了能在R中存放如此大的数据，IRanges使用一个小技巧来压缩这些序列。例如444333222就可以表示为3个4,3个3，3个2。这个压缩方法称为*run-length encoding, RLE* 。并且压缩后的数据支持原来的所有常规的R向量操作，如截取、算数、summary,和数学函数

```r
x <- as.integer(c(4, 4, 4, 3, 3, 2, 1, 1, 1, 1, 1, 0, 0, 0,0, 0, 0, 0, 1, 1, 1, 4, 4, 4, 4, 4, 4,4))
xrle <- Rle(x) # 压缩为run coding length

integer-Rle of length 28 with 7 runs
  Lengths: 3 2 1 5 7 3 7
  Values : 4 3 2 1 0 1 4

as.vector(xrle) #解压缩
runLength(xrle) # 仅获取lengths
runValue(xrle) # 仅获取values
```

我们可以在converage()中遇到run-length encoded values. converage()输入一组ranges, 返回以Rle对象表示的覆盖度.

```r
set.seed(0)
rngs <- IRanges(start =sample(seq_len(60), 10), width=7)
names(rngs)[9] <- 'A'
rngs_cov <- coverage(rngs)

rngs_cov > 3 # 找到覆盖度大于3的区域
rngs_cov[as.vector(rngs_cov) >3] # 提取覆盖度大于3的值
rngs_cov[rngs['A']] #标记为A的区域的覆盖度
mean(rngs_cov[rngs['A'])
```

在基因组学中有大量的分析工作，例如对一组ranges（重复区域，蛋白编码区域，低重组区等）计算一些序列的描述性统计量（coverage, GC含量，核酸多样性等）。这些操作在我们使用IRanges中的methods会选择特别琐碎，不过可以用GenomicRanges中相似方法来处理。

## 使用slice()从RLE序列创建ranges

define new ranges: taking coverage and defining ranges correspnding to extremely high-coverage peaksm or a map of per-base pair recombination estimates and  defining a recombinationa hotspot region.
例如找出覆盖率低于2的ranges.

```r
min_cov2 <- slice(rngs_cov, lower=2)
min_cov2
Views on a 63-length Rle subject

views:
    start end width
[1]    16  18     3 [2 2 2]
[2]    22  22     1 [2]
[3]    35  39     5 [2 2 2 2 2]
[4]    51  62    12 [2 2 2 3 3 3 4 3 3 3 2 2]

rangs(min_cov2) # 提取views中的range
```

## IRanges进阶：Views

Views通过整合sequence vector和ranges,使得对根据特定ranges的sequence vector的聚合操作更加简单。这与之前所提到分组描述性统计分析类似，只不过这里的组别为ranges。
例如，我们可以对之前使用slice()创建的views进行描述性统计分析，使用的函数为viewMeans(), viewMaxs(),和viewApply(可以使用任意函数）

```r
viewMeans(min_cov2)
viewMaxs(min_cov2)
viewApply(min_cov2, median)
```

还有其他类似的函数，见help(viewMeans)。
同样还可以根据window/bin对序列进行描述性分析。首先是在序列中创建一组windows ranges,然后在分析.

```r
length(rngs_cov)
bwidth <- 5L
end <- bwidth*floor(length(rngs_cov / bwidth))
windows <- IRanges(start= seq(1,end, bwidth), width = bwidth)
head(windows)
cov_by_wnd <- View(rngs_cov, windows)
head(cov_by_wnd)
viewMeans(cov_by_wnd)
```

命令行有一个工具叫做bedtools，能够实现如上代码所做的事情。