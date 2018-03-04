# 如何分析遗传数据

通过高通量测序拿到了大量样本的SNP信息后，下一步该怎么办？不妨看看这张<http://www.molecularecologist.com/>在2018年1月份的一篇文章《Molecular ecology, the flowchart》里的套路图吧

![分析套路](http://oex750gzt.bkt.clouddn.com/18-2-19/37890748.jpg)

## 描述群体多样性

任何一个遗传群体都是由它所包含的各种基因型所组成的，在一个群体内某特定基因型所占的比例，就是基因型频率(genotype frequency)。基因型是每代在受精过程中由父母所具有的基因所组成，**它是描述群体遗传结构的重要参数**。在一群体内某特定基因座某一特定等位基因占该基因座等位基因总数的比率，即为等位基因频率(allel frequency),或称基因频率。**等位基因频率是决定一个群体遗传特性的基本因素**。当环境或遗传结构不变时，等位基因频率不会改变。

哈迪-温伯格(Hardy-Weinberg)定律:在一个完全随机交配的群体内，如果没有其他因素（如突变、选择、迁移、遗传漂变等）干扰时，则等为基因频率及3中基因型频率始终保持一致，各代不变。

核酸多样性（Nucleotide diversity）是分子遗传学中衡量群体多态性的一个概念，最早由Li和Nei在1979年提。计算公示如下

![核酸多样性](https://wikimedia.org/api/rest_v1/media/math/render/svg/be2956df9d2756a4f051f2516938d4831fcd3771)

其中xi和xj表示第i和第j条DNA序列的频率，pai\_ij则是第i和第j条序列每个核苷酸位置上的差异。

## 检查群体结构

### Fst

Fst(Fixation index,固定指数)是衡量因遗传结构(genetics structure)而引起群体差异(population differentiation), 是Sewall Wright提出的F-统计值的特例，可以说是目前群体遗传学最常用的统计值。

> Sewall Wright是现代群体遗传学的奠基者之一，其余两人是Ronald Fisher（费希尔）和 J.B.S. Haldane。群体遗传学好多概念就是他提出的，比如说近交系数，遗传漂变等。



[Fst和Gst](http://www.molecularecologist.com/2011/03/should-i-use-fst-gst-or-d-2/)