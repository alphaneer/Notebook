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

对于高通量测序全（简化）基因组测序而言，目前的公司测序体系比较成熟，提供的原始数据质量都比较高，不需要过多处理。目前150bp PE read只需要去掉低质量序列和接头即可，对应软件为trimmomatic和TrimGalore。

```bash


```

着丝粒附近假阳性SNP位点高，且Q-values低