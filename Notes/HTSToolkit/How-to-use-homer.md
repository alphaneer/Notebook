# 如何使用HOMER找Peak

找Peak目前的主流软件，一个是MACS2，另一个就是这篇介绍的HOMER。HOMER最早是用来寻找motif, 后来是整合了ChIP-seq, GRO-seq, RNA-seq, DNase-seq, Hi-C等数据分析的工具，功能更加强大，也一直有在维护更新。软件对应的参考文献是"Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities"。

目前HOMER支持大多数主流的模式物种，同时允许提供FASTA格式的参考基因组和GTF格式的注释文件自行构建

## 软件安装

官方提供了`configureHomer.pl`用于安装和升级HOMER

```bash
cd ~/opt/biosoft
mkdir HOMER && cd HOMER
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install homer
```

如果要做差异peak分析，你需要额外安装两个R包,DESeq2和edgeR。比较简单的方法是通过bioconda

```bash
conda install wget samtools r-essentials bioconductor-deseq2 bioconductor-edger
```

此外通过`perl configureHomer.pl -list`了解可供下载的数据，选择目标物种下载。

```bash
perl configureHomer.pl -install arabidopsis-p
perl configureHomer.pl -install arabidopsis-o
perl configureHomer.pl -install tair10
```

## 用findPeaks找peak

HOMER使用findPeaks进行peak calling, 命令的常规使用方法为

```bash
findPeak <实验组文件路径> -style <运行模式> -o <输出文件名> -i <对照组文件路径>
```

findPeaks提供了7种运行模式，如下所示

- factor: 适用于单个结合位点的ChIP-seq或DNase-seq, 用于预测转录因子的结合位点，蛋白和DNA的结合位置。
- histone: 不同的组蛋白标记的ChIP-seq试验
- super: 超级增强子模式
- groseq: 适用于链特异性的GRO-seq
- tss: 从5'RNA-seq 或 5'GRO-seq数据中寻找启动子
- dnase: 适用于DNas-seq peak calling
- mC: 适用于DNA甲基化数据分析