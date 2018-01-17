---
title: 利用ggbio对生物数据进行可视化
tags: R, bioconductor, visualization
notebook: 工具笔记
---
# 利用ggbio对生物数据进行可视化

ggbio是图形语法(graphic grammar)在生物数据领域上的拓展，任何ggbio得到的结果都能与ggplot2进行互动，而不是简单封装。

## tracks

ggbio对ggplot2的一个拓展是其提供了`tracks`用于以不同轨的方式比较数据，比如说比较两个不同的时间序列。

```r
library(ggbio)
df1 <- data.frame(time = 1:100, score = sin((1:100)/20) * 10)
p1 <- qplot(data = df1, x = time, y = score, geom = "line")
df2 <- data.frame(time = 30:120, score = sin((30:120)/20) * 10, value = rnorm(120 -
30 + 1))
p2 <- ggplot(data = df2, aes(x = time, y = score)) + geom_line() + geom_point(size = 4,
aes(color = value))
lst <- list(time1 = p1, time2 = p2)
tracks(lst)
```

![时间序列轨](http://oex750gzt.bkt.clouddn.com/18-1-17/60903679.jpg)

tracks支持多种方式定制

- height: 定义某个track的高度
- bcColor: 定义某个轨道背景颜色
- labeled: 对某个轨道的标签进行命名
- fixed: 控制标度的固定还是不固定
- mutable: 控制作图结果能否使用`+`更改tracks

此外tracks还拓展了ggplot2的主题，以`theme_tracks_*`方式命名.

最后这些修改方式都可以通过`reset`一键还原到之前`backup`的地方。

## mold

ggbio支持Bioconductor的几个核心数据结构的可视化, 如IRange和GenomicRange, 而ggplot2则是主要处理data.frame.为了方便衔接，作者开发了`mold`将Biconductor的核心数据结构转换成data.frame. `mold`相对于`as.data.frame`这种简单粗暴的转换，保留了原来数据结构的更多信息，创建了更多变量统计值。目前支持eSet, GRanges, IRanges, GRangesList, Seqinfo, matrix, Views, ExpressionSet SummarizedExperiment, Rle, RleList。

PS: `mold`函数在`biovizBase`包中。

## ggbio拓展的图形语法

判断一个人是否懂ggplot2，只要问他"什么是图形语法"即可

> 一张统计图形就是从数据到**几何对象**(geometric object, 点、线、大小等)的**图形属性**(aesthetic attribute, 颜色、形状、大小)的一个映射。此外，图形中还可能包括数据的**统计变换**(statistical transformation), 最后绘制在某个特性的**坐标轴**(coordinate system)中，而**分面**(facet)则可以用来生成数据不同子集的图形。

因此真正理解图形语法的人，看了下面这张图就会用ggbio了。

![](http://oex750gzt.bkt.clouddn.com/18-1-17/79879573.jpg)

有些细节不懂，就可以翻翻<http://www.bioconductor.org/packages/release/bioc/manuals/ggbio/man/ggbio.pdf>

## 几个案例

### 基因结构图

这部分操作涉及到`AnnotationHub`和`GenomicFeatures`两个R包，这些在我的[用Bioconductor对基因组注释](https://www.jianshu.com/p/ae94178918bc)有过比较详细的介绍。如果理解了那两个R包的使用，那么起始用起来也就很简单，无非就是构建一个ggbio能用的GRanges对象而已，步骤如下：

1. 使用select根据基因ID找到CDS ID
1. 使用cds根据CDS ID筛选区域
1. 对区域使用`geom_arrowrect`作图

```r
# 加载拟南芥TxDb
library(AnnotationHub)
ah <- AnnotationHub()
TAIR_tx <- ah[['AH52247']]

# 1. 根据基因ID找到CDS ID
keys <- c('AT3G27920')
cds_region <- select(TAIR_tx, keys =keys, columns = 'CDSID', keytype = 'GENEID')
# 2. 根据CDS ID创建GenomicRanges
cds_gr <- cds(TAIR_tx, filter = list(cds_id = cds_region$CDSID))
# 3. 作图
library(ggplot2)
ggplot(cds_gr) + geom_arrowrect(fill='black') + theme_clear() + theme_alignment()
```

![](http://oex750gzt.bkt.clouddn.com/18-1-17/49948970.jpg)
