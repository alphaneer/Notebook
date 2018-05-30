# GATK 4.x版本常用工具介绍

## 增加或替换read group

如果一个BAM文件中没有RG这一项，那么GATK相关程序运行时就会出现问题。报错还好，你还能检查，如果没有报错，一直等到软件运行结束后你才发现结果异常，就会浪费时间。

GATK4.x版整合了PICARD的`AddOrReplaceReadGroups`,使用方法如下

```bash
~/opt/biosoft/gatk-4.0.3.0/gatk AddOrReplaceReadGroups -I input.bam \
    --RGLB 文库编号  \
    --RGPL 测序平台 \
    --RGPU 平台单元 \
    --RGSM 样本名字 \
    -O output.bam
```

GATK考虑到同一个样本测序量很大，可能多次建库，用不同的测序平台，然后每个样本测序还会用条形码(barcode)进行区分。

```bash
~/opt/biosoft/gatk-4.0.3.0/gatk AddOrReplaceReadGroups \
    -I input.bam \
    --RGLB lib1  \
    --RGPL illumina\
    --RGPU unit1 \
    --RGSM sampleA \
    -O output.bam
```

建议写比对参数的时候，不要偷懒，把RG加上, 比如说Bwa和bowtie2

```bash
bowtie2 --rg-id "文库编号" --rg "SM:样本名" --rg "PL:ILLUMINA"
bwa mem  -R '@RG\tID:文库编号\tSM:样本名\tPL:ILLUMINA'
```

## 拆分VCF文件

SplitVcfs (Picard)