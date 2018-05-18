# 如何分析芯片数据

我最早接触的高通量数据就是RNA-seq，后来接触的也基本是高通量测序结果而不是芯片数据，因此我从来没有分析过一次芯片数据，而最近有一个学员在看生信技能树在腾讯课堂发布的课程[GEO数据库表达芯片处理之R语言流程](https://ke.qq.com/course/286407?tuin=1ae7bc83)遇到了问题问我请教，为了解决这个问题，我花了一个晚上时间学习这方面的分析。 **注**:这篇文章不会介绍R语言的安装和使用，也不会介绍GEO数据库的构造

## 数据的获取

数据获取有两种方式，R包**GEOquery**解析和手动下载。其中前面一种最方便，完成了手动数据下载和Bioconductor常见数据结构`ExpressionSet`的构造，关于这个数据结构的具体介绍看Bioconductor的介绍或者视频，简言之，就是用于存放 **实验信息**, **分组信息** 和 **表达信息**, 方便后续调用。

```bash
library(GEOquery)
gset <- getGEO("GSE13535", GSEMatrix =TRUE, AnnotGPL=TRUE )
show(gset)
```

![ExpressionSet](http://oex750gzt.bkt.clouddn.com/18-5-16/48674991.jpg)

一般而言GEOquery解析都是首选，除非你提供的GSE号还没被GEOquery记录或者说网络速度感人，以及你不觉得别人提供的矩阵是你所需要的，你才会决定去手工下载。分为两种情况，一种是下载赛默飞的下机原始数据格式CEL，一种是下载单个样本表达量向量或者含有所有样本的表达量矩阵。

![数据下载](http://oex750gzt.bkt.clouddn.com/18-5-16/95577053.jpg)

先说第一种，可以直接点击http下载到tar打包的数据, 然后解压缩得到所有的CEL文件

```r
setwd("F:/Project/GEO_project/")
library(affy)
affy.data <- ReadAffy()
length(affy.data)
# 13
eset.rma <- rma(affy.data)
exprSet <- exprs(eset.rma)
write.table(exprSet, "expr_rma_matrix.txt", quote=F, sep="\t")
```

- ReadAffy: 读取当前文件下的CEL格式文件，同时第一次还会从bioconductor上下载hugene10stv1用来注释cel文件。
- rma: 基于robust multi-arrary average(RMA)算法衡量表达量，从而将AffyBatch对象转换成**ExpressionSet**
- exprs: 获取`ExpressionSet`中的表达量矩阵
- write.table: 将表达量矩阵信息保存到本地

然后是第二种，以所有样本的表达矩阵为例，可以用浏览器到<ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42589/matrix/>下载，如果你会用Linux的话，可以用`wget -4 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42589/matrix/GSE42589_series_matrix.txt.gz`, 才1.7M。解压缩这个文件后，有一个txt文件, 这个txt分为两个部分。第一个部分是以"!"开头的样本的所有信息，如实验平台、处理、以及分组等信息。第二个部分则是后面的表达量信息，

![Series Matrix Files](http://oex750gzt.bkt.clouddn.com/18-5-16/54325973.jpg)

```bash
expr.df <- read.table(file = "GSE42589_series_matrix.txt", header =TRUE,
                      comment.char = "!", row.names=1)
```

可以从这个角度理解这三种方法： 最开始得到的都是CEL文件，CEL文件需要一系列的步骤才能转换成表达矩阵，例如去除批次效应、质控和过滤等，得到的表达矩阵在上传时会增加元数据信息（处理方法、分组信息），就成为我们下载的`GSEXXXX_series_matrix.txt.gz`. 通过手工解析加R语言简单操作得到了R语言中的数据框(data.frame)， 而GEOquery能够帮助我们完成下载和解析这两个步骤。

三者的优先级为：GEOquery > 手工下载表达量矩阵文件 > 手工下载原始的CEL文件。

## 使用limma进行差异表达分析

limma的核心函数是lmFit和eBayes， 前者是用于线性拟合，后者根据前者的拟合结果进行统计推断。

lmFit至少需要两个输入，一个是表达矩阵，一个是分组对象。

**表达矩阵**必须是matrix类数据结构，每一列都是存放一个样本，每一行是一个探针信息或者是注释后的基因名。这里就是向我提问的人出错的原因，他在读入数据时，read.table少了参数，`row.names= 1`，导致第一列是探针信息。

```r
# 使用GEOquery
exprSet <- exprs(gset[[1]])
# 基于matrix
expr.df <- read.table(file = "GSE42589_series_matrix.txt", header =TRUE,
                      comment.char = "!", row.names=1)
# 从cel文件开始
exprSet <- exprs(eset.rma)
```

**试验设计矩阵**: 没有试验设计矩阵对象，limma就不知道如何比较。分组数据可以手工从之前的matrix.gz整理，整理到一个excel，然后用R读取，或者就是直接从Geoquery的结果中解析。

```r
pData <- pData(gset[[1]])
view(pData)
```

![GEOquery解析的信息](http://oex750gzt.bkt.clouddn.com/18-5-16/94911404.jpg)

其中title部分告诉了我们分组信息，2小时和18小时，每个时间段又有vehicle control, PE1.3 embolized, PE2.0 embolized，也就是2x2的双因素试验设计, 我们可以现在R语言里构建实验设计的数据框。

```r
sample <- pData$geo_accession
treat_time <- rep(c("2h","18h"),each=11)
treat_type <- rep(rep(c("vehicle_control","PE1.3_embolized","PE2.0_embolized"), c(3,4,4)),
                  times=2)
design_df <- data.frame(sample, treat_time, treat_type)
```

根据Limma的使用手册的"9.5 Interaction Models: 2 X 2 Factorial Design"进行手续的分析。这里仅仅展示单个因素的分析过程，多个因素看文档依葫芦画瓢就行。

构建单因素试验设计矩阵，进行线性拟合

```bash
TS <- paste(design_df$treat_time, design_df$treat_type, sep=".")
TS
TS <- factor(TS, levels = unique(TS))
design <- model.matrix(~0+TS)
fit <- lmFit(exprSet, design)
```

然后根据我们要回答的问题，来设置比较对象。比如不同时间段下控制组哪些基因发生了差异报答，处理18小时后，处理组相对于对照组有哪些基因发生差异表达，也就是做多次t检验。

```r
cont.matrix <- makeContrasts(
  vs1  = TS18.vehicle_control-TS2.vehicle_control, # 对照组在前后的差异表达基因
  vs2  = TS18.PE2.0_embolized-TS2.PE2.0_embolized, # PE2.0处理前后的差异基因
  vs3  = TS18.PE1.3_embolized-TS2.PE1.3_embolized, # PE1.3在处理前后差异基因
  # 处理18小时候，PE2.0相对于对照变化的基因再与PE1.3与对照的差异比较
  diff = (TS18.PE2.0_embolized-TS18.vehicle_control)-(TS18.PE1.3_embolized-TS18.vehicle_control),
  levels = design
)

fit2 <- contrasts.fit(fit, cont.matrix)
results <- decideTests(fit2)
```

最后的结果可以用韦恩图展示`vennDiagram(results)`

## 更多分析

找到的差异表达基因后续要做GO/KEGG分析，可以在生信技能树公众号中搜索，或者阅读原文购买视频。