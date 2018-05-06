# 如何过滤freebayes得到原始vcf文件

变异检测算法的核心就是从尽可能找到真实的变异，降低假阳性。尽管目前测序仪器的准确性可以达到99.999%，似乎很高的样子，但是对于高通量测序而言，这意味着在100,000个碱基中就可能出现一个错误，那么freebayes如何保证自己结果的可靠性？

freebayes基于贝叶斯公式

![贝叶斯公式](http://oex750gzt.bkt.clouddn.com/18-4-25/95369344.jpg)

简单的说，当一个变异如果只出现在一条链上，或者是某一个位置上，那么这个位点很有可能是高通量测序时引入的偏误。

![可能位点](http://oex750gzt.bkt.clouddn.com/18-4-25/258200.jpg)

先验模型并不能解决所有错误，freebayes初步会得到海量的变异位点，这肯定是不能直接用于最后分析，需要进一步过滤。过滤有两种策略，一种是硬过滤(hard filter)，一种则是使用机器学习的方法，比如说支持向量机。

Hard filters的策略很简单，就是按照**我们所认为的好**去过滤，

- 这个变异的信度要高, QUAL>N
- 有足够多的深度支持，DP >N
- 变异应该出现在两条链上, SAF >0 & SAR>0
- 变异出现在read的中部, RPL>0 & RPR >0

> RPL(Reads Placed Left), RPR(Reads Placed Right)
> SAF(Number of alternate observations on the forward strand), SAR(Number of alternate observations on the reverse strand)

由于每一个物种基因组性质都不太一样，那么应该设置什么样标准比较好呢？

如果是自然变异, 那么在大多数生物中，转换(transitions, ts, A-T<->G-C)的发生的概率颠换应该大于颠换(transversion, tv,T-A<->G-C). 在人类中,ts/tv约等于2, 在线粒体中, ts/tv有可能大于20。 tv/ts信息可以用`vt peek`看。

使用机器学习的方法有点难度，需要你提供一个高信度的变异集用来训练分类器，属于比较高级的模块。

参考资料

- <https://github.com/ekg/freebayes>
- <http://ddocent.com/filtering/>
- <https://github.com/ekg/alignment-and-variant-calling-tutorial>