# BSA流程学习的二八法则

## 实战部分

数据变换

fastq $\to$ SAM $\to$ BAM $\to$ **VCF** $\to$ Results

fastq $\to$ SAM : BWA, Bowtie2, SOAP

SAM $\to$ BAM : SAMTools

BAM $\to$ **VCF** : GATK,  SAMTools + VCFTools, freebayes

BWA $\rightarrow$ SAMTools $\to$ VCFTools $\to$ R/SHOREmap/MutMap



```shell
#!/bin/bash


```





