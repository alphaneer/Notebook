---
title: 全基因组重测序变异检测
tags: 框架, framework
notebook: 分析流程
---

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

# 全基因组重测序变异检测

## 项目结构

项目结构结构分为三个部分，运行脚本，原始数据和分析结果

```bash
project-snp
|-- README.md
|-- run
|-- data
|   |-- raw_data
|   |-- clean_data
|-- analysis
    |-- alignment
    |-- mpileup
    |-- variant
        |-- snp_indel
```

## 分析流程和对应代码

### 原始数据预处理

对于高通量测序全（简化）基因组测序而言，目前的公司测序体系比较成熟，提供的原始数据质量都比较高，不需要过多处理。目前150bp PE read只需要去掉低质量序列和接头即可，对应软件为trimmomatic(Java)

```bash


```

着丝粒附近假阳性SNP位点高，且Q-values低

### 变异注释

GATK在变异检测的时候会自动计算一些注释信息，常用注释模块解释如下

- BaseQualityRankSumTest(BaseQRankSum): 通过秩和检验，分析支持参考位点和支持变异位点的碱基质量是否存在差异性。最好是0，如果是负数则表示支持参考位点的read碱基质量比较高，反义亦然。
- RMSMappingQuality(MQ)：估计了支持变异位点的所有reads的比对质量，60封顶。
- ReadPosRankSumTest(ReadPosRankSum): 通过秩和检验，检验支持变异位点的碱基在各自的read上是不是存在位置偏好性。毕竟如果支持变异的碱基都在reads的末端，有很大可能是假阳性。理想状态是0，如果是负数，表明变异的碱基更多位于末端。
- FisherStrand(FS): 判断是否存在链偏好性，越小越好。

得到的VCF文件还能进一步用`VariantAnnotator`注释，案例代码

```bash
java -Xmx4g -jar GenomeAnalysisTK.jar -T VariantAnnotator \
    -R referecne.fa -I input.vcf -L input.vcf \
    -A Coverage # 注释深度
    -o output.vc
```

如果想知道一个碱基的变异效应，可用snpEff注释后，使用`VariantAnnotator`进行整合

```bash
# annotation with snpEff 
java -Xmx4g -jar snpEff.jar \
  -v \
  -o gatk \
  referece_db \
  input.vcf \
  > input.ann.gatk.vcf
# integrated with GATK VariantAnnotator
java -Xmx4g -jar $HOME/tools/gatk/GenomeAnalysisTK.jar \
    -T VariantAnnotator -R reference.fa
    -A SnpEff \
    --variant input.vcf \
    --snpEffFile input.ann.gatk.vcf \
    -L input.vcf \
    -o output.gatk.vcf
```

> 参数-L可以选择注释的范围，例如 -L chr1:100-200 表明对chr1的100~200范围进行注释。

### 变异过滤​

变异过滤可以选择的工具比较多，例如vcftools, plink甚至自己写脚本都是可行。

假如你需要按照如下规则进行过滤：

- quality(Q) > 30
- mapping quality(MQ) >20
- quality-by-depth-ratio (QD) > 10
- ReadPosRankSum > -8.0
- depth coverage(DP) > 3
- probability of strand bias(FS) < 10(snp), 对于indel则是200
- 10bp 内不能多于三个snp