---
title: 如何对得到的变异进行过滤
author: xu zhougeng
tags: pipeline
notebook: 分析流程
---
# 如何对得到的变异进行过滤

变异检测（variant calling）这一步得到的VCF/BCF文件中有着大量的变异位点。由于工具本身的算法特性，同一批数据里得到的变异结果有着不同程度上假阳性，因此需要先进行过滤，才能用于后续的分析。过滤分为两种，一种是GATK的VQSR，通过GMM模型训练分类器进行区分，另一种是“hard filter”，直接根据INFO里的Tag对结果进行过滤。前者需要比较大对群体，需要较为完善的数据库，对于一般的动植物，可能没有那么多群体和完备的数据库用来训练分类模型，因此更多采用“hard filter”的方式。

对于“hard filter”，也就是硬过滤，也有两种类型过滤方法，一种是在VCF的header部分添加FILTER信息，FILTER列则以PASS和对应的标准表示，之后还要根据是否PASS来保留相应的位点。另一种则是一步过滤。

变异过滤可以用到的工具为：

- bcftools
- vcftools
- SnpSift
- Picard

我们主要需要了解应该基于什么样的标准进行过滤，为了实现这种方法的过滤，需要用到上面的哪一个工具。

## GATK推荐标准

GATK的HaplotypeCaller找到的变异可以按照GATK推荐的标准进行过滤，即对于SNP为`"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"`，对于InDel为`"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"`,用的是JEXL表达式，里面的参数含义为

- QualByDepth (QD)： 变异质量QUAL的另一种表示方法，越大越好
- FisherStrand (FS)： 越接近0越好，正常的变异应该均衡在两条链上检测出
- RMSMappingQuality (MQ)： 所有样本比对质量的均方根，越大表示比对质量高
- MappingQualityRankSumTest (MQRankSum): 越接近0越好，变异位点和非变异位点
- ReadPosRankSumTest (ReadPosRankSum)：越接近0越好，用来剔除那些在特异性出现在read末端的偏差
- StrandOddsRatio (SOR)：越低越好，和FS类似

在GATK4版本中，picard已经整合到gatk中，假定得到的所有变异文件命名为raw\_variants.vcf,参考基因组命名为referece.fa,snp的过滤代码如下：

```bash
# 筛选snp
~/biosoft/gatk-4.0.1.2/gatk SelectVariants \
    -R reference.fa \
    -V raw_variants.vcf \
    -select-type SNP \
    -O raw_snps.vcf
# snp 过滤
~/biosoft/gatk-4.0.1.2/gatk VariantFiltration \
    -R reference.fa \
    -V raw_snps.vcf \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "gatk_snp_filter" \
    -O filtered_snps.vcf
# 可以用grep简单了解下多少snp会被过滤
grep -c 'gatk_snp_filter' filtered_snps.vcf
```

对应的InDel过滤代码如下

```bash
# 筛选indel
~/biosoft/gatk-4.0.1.2/gatk SelectVariants \
    -R reference.fa \
    -V raw_variants.vcf \
    -select-type INDEL \
    -O raw_indels.vcf
# snp 过滤
~/biosoft/gatk-4.0.1.2/gatk VariantFiltration \
    -R reference.fa \
    -V raw_snps.vcf \
    --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
    --filter-name "gatk_indel_filter" \
    -O filtered_indels.vcf
# 可以用grep简单了解下多少snp会被过滤
grep -c 'gatk_snp_filter' filtered_indels.vcf
```

## SNP过滤参数

对于动植物研究而言，最常用的分子标记就是SNP（SNP在当前语境指的是和参考基因组不同的单的碱基差异），因此这里说点和SNP过滤的参数。

**SNP间距**：一般来说，自然界中基因突变率比较低，并且发生的位置相对随机。这意味着对于同一个物种，在一个比较小的区域，比如说1000 bp里出现了10个以上的snp，就有很大可能是假阳性。处理这类位点，可以自己写脚本，但是GATK4/VariantFiltration提供了`--cluster-size,-cluster:Integer`和`--cluster-window-size,-window:integer`。

**InDel附近的SNP**：Indel附近的snp通常不太准确，不过这个问题HaplotypeCaller已经考虑到了，它会排除离InDel15bp以内的SNP。如果你需要进一步剔除50bp内或者稍微再远一点的SNP，可以使用`bcftools filter --snpGap`.

**SNP深度**： 对于全基因组深度测序而言，我们根据理论的测序深度剔除过高或者过低覆盖的碱基。过低的位点从逻辑上比较好理解，过高的位点可能是由于参考基因组里面本来没有的拷贝数变异和旁系同源序列引起。

## 再说几句

过滤其实就是取舍的过程，降低假阳性的同时也会提高假阴性。当然在高通量测序中，我们的烦恼可能更多是来自于噪声太大，也就是位点太多。

## 参考文献

- <http://kaopubear.top/2018-01-31-callvariantfilter.html>
- <https://software.broadinstitute.org/gatk/documentation/article.php?id=1255>
- <https://software.broadinstitute.org/gatk/documentation/article.php?id=2806>
- <https://link.springer.com/article/10.1186/s13073-014-0089-z>