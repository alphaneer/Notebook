# 高通量测序数据得到的SNP如何构建高密度遗传图谱

## 作图群体

- 暂时群体: F2, BC
- 永久性群体: RIL, CSSL, DH
- 自然群体

## SNP标记

高通量SNP分子标记的可能问题, 假阳性高

测序问题：

- 碱基错误
- 深度过高或者过低

选择测序深度在平均深度1/3到2倍的标记

基因型问题:

- 数据缺失比例： 亲本多态的位点在后代中未必存在
- 过多的单交换和非预期双交换(unexpected double recombinants)

![双交换VS单交换](http://oex750gzt.bkt.clouddn.com/18-4-27/57197203.jpg)

- 等位基因变换(alleles switching)
- 杂合度,也就是偏分离程度(segregation distortion)：在连锁图谱中，偏分离分子标记往往扎堆形成偏分离区(segregation distorted regions), 通常过滤p<0.001的点。

以上这些问题都会导致遗传图谱的大小预测不准确，分子标记位置不对。综上，高通量测序得到的SNP标记还需要得到强有力的过滤和筛选, 比如说下面这几篇文章

- Waseem et.al 2017年的文章发现了6120个分子标记，其中1303个标记明显的偏向一方亲本，并研究这些标记在基因组上偏分离的模式。
- Zhang et.al 2016年找到了双亲中40,372个多态标记, 但是只用了2540(1756非冗余)标记做了遗传图谱, 平均距离0.88 cM.
- Yang et.al 2017找到了双亲中21,584个多态标记，最后只使用了3655个标记, 0.81cM.

也就是说初步过滤后的亲本间多态标记最后也只会有不到20%的比例用于遗传作图。并且从这些文章中，可以总结出如下过滤规律:

- 分子标记在后代的缺失比例不能过高，要低于25%
- 用于构图的分子标记的偏分离不能过于严重, 卡方检验低于0.001
- 分子标记的相似度要低于100%

尽管看起来分子标记很多，但其实SNP有很多是共分离的标记，可以进行合并。在群体创建时，交换数目就已经确定了，F2的交换少，RIL多。

## 构图软件

算法分类：

- multipoint-likelihood maximizaition: MapMaker, CRI-MAP, CarthaGene, R/qtl
- two-point statistics: GMendel, JoinMap, RECORD

区间作图在理论上更加有优势但是速度上慢

分群

- Joinmap(最权威)
- lcimapping
- mstmp, 适合高通量, 最小生成树法, 适用于DH, BC1, Hap, RIL(n)
- highmap, 打包mstmap和矫正程序

蚁巢算法

分群后分别做遗传图

## 参考资料

- <https://biocyclopedia.com/index/genetics/linkage_and_crossing_over_in_diploid_organisms_higher_eukaryotes/coupling_and_repulsion_hypothesis.php>
- Genetic Mapping in the Presence of Genotyping Errors
- Genotyping-by-Sequencing Derived High-Density Linkage Map and its Application to QTL Mapping of Flag Leaf Traits in Bread Wheat
- Development of a high-density linkage map and mapping of the three-pistil gene (Pis1) in wheat using GBS markers
