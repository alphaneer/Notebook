# 使用二代测序寻找T-DNA插入

## 通过数据模拟和可视化方式预探索

为了解基因组存在T-DNA插入时，即基因组构成为AC而样本基因组为ABC的情况得到的测序结果在序列比对的时候的可能情况，因此需要先要使用模拟数据进行探索。

第一步：截取拟南芥参考基因组

## 实际流程

第一步：将T-DNA序列附加到参考基因组，并建立索引

```bash
cat t-dna.fa >> reference.fa
bwa index reference.fa
```

第二步：将PE序列比对到附加了T-DNA序列的参考序列

```bwa
bwa mem -t 10 reference.fa read_1.fq read_2.fa | samtools sort > aligned.bam
samtools index aligned.bam
```

第三步：了解一下序列比对情况.结果分为四列:染色体名，染色体长度，比对序列read数，未比对序列数

```bash
samtools idxstats aligned.bam
```

第四步：从比对序列中提取出如下几类read

- PE reads其中一条比对到T-DNA，另一条比对到reference
- PE reads尽管全部比对到T-DNA序列，但是至少有一条read存在soft clip和hard clip现象
