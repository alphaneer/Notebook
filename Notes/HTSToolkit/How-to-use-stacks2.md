# 如何使用stacks2分析RAD-seq

Stacks2.0(后面无特殊说明简称stacks)是RAD-seq数据的分析流程，支持有参考和无参考基因组的RAD-seq分析。目前支持如下类型的数据

- 文库类型: 原始RAD-seq, ddRAD-seq, 2bRAD, GBS, CRoPS
- 测序类型: illumina,  Ion Torrent
- 短读类型: 双端测序(paired-end), 单端测序(single-end). 对于原始的RAD-seq, stacks会尝试进行组装， 对于ddRAD-seq或GBS这类，双端会形成两个loci。最后这两个loci合并成一个loci, 如果这个loci不重合，中间用N补全。

stacks工具分为三个部分:

- 原始数据处理: process\_radtags, process\_shorreads, clone\_filter, kmer\_filter
- 核心部分: ustacks, cstacks, sstacks, tsv2bam, gstacks, populations
- 执行脚本: denovo\_map.pl, ref\_map.pl

如果公司提供的是去掉barcode数据(clean data)，那么原始数据处理这个部分可以省去。或者使用其他工具处理，但是要注意以下几点，不然后面会报错。

- 双端测序的文件名一定要是: xxx.1.fq.gz, xxx.2.fq.gz
- 双端测序的read ID一定要是: xxxx\_1, xxxx\_2 或 xxxx/1, xxxx/2

stacks的核心在于他开发的6个工具，功能如下

- ustacks: 将输入的短读序列比对到准确匹配的stacks，然后比较这些stacks通过最大似然法识别出潜在的座位(loci)和SNP
- cstacks: 整合ustacks找到的各个样本的位点，构建出一致性座位(loci), 最后合并等位基因形成catalog(等位基因一览表)。在遗传杂交群体中，则会根据杂交的亲本构造出catalog, 来得到一组杂交后代所有可能的等位基因
- sstacks: ustacks得到的可能座位(loci)进一步在catalog(来自于cstacks)搜索。如果是遗传群体, 会用后代的stacks匹配catalog来确定那个后代包含哪个亲本的等位基因。如果是常规群体，则是使用所有群体。
- tsv2bam: 该程序用来数据格式转换，原本的数据以样本主导，转换后就是以座位(locus)导向。并且当PE数据可用时，这个程序还会将每个单端座位(loci)相关的双端数据进行从头组装
- gstacks: 这个工具有两个模式，无参分析时会用双端测序组装成contig然后合并成单个座位(locus)，接着将read比对到座位上寻找SNP。 有参分析，则是直接使用比对并排序好的数据构建座位。
- populations: 根据以上基本得到的数据进行群体遗传学分析。

如果样本量不大，或者是单节点计算，官方提供了两个脚本整合了上面几个工具，分别是**denovo\_map.pl**和**ref\_map.pl**, 其中**denovo\_map.pl**在对参数进行调优时特别实用, 毕竟只会用到20个样本。

## denovo\_map.pl使用说明

`denovo_map.pl`运行时需要准备

```bash
denovo_map.pl --sample
```

## 参数选择

Stacks主要需要调整的参数为`ustacks`的`-M, -m`和`cstacks`的`-n`, 其中`ustacks`的`M`控制两个不同样本等位基因（allele）之间的错配数，`m`控制最少需要几个相同的碱基来形成一个堆叠（stack).。而cstacks的`n`和`ustacks`的`M`等价。最后一个比较复杂的参数是能否允许gap(--gap)

这三个参数的设置取决于你RAD-seq的几个主要特征：

- 生物学背景如倍性、多态性和需要检验的生物学假设
- 不同RAD(RAD, ddRAD, 2bRAD)和所使用的限制性内切酶导致的原始短读长度和数量，酶切位点、覆盖度，以及测序平台本身的系统误差
- 文库质量: DNA降解程度, 外界污染

在"Lost in parameter space: A road map for Stacks"一文中给了三个参数的选择建议

| 参数 | 范围 | 决策1 | 决策2 | 决策3 | 其他建议|
| ----|  --- | ----- | ----- | -----| -------|
| m  默认3 | 3 ~7 | 覆盖率 < 15时提高 | 有污染离子时提高 | 系统发育研究时提高 | 当m>6时不要使用第二个reads |
| M  默认2 | 1 ~8 | 多态性高时提高    | 基因组分离程度高时提高 | 重复序列或多倍体降低 | 过高的M会导致`popluations`过滤旁系同源位点。在读长大于250时需要重新调整参数|
| n 默认1 | =M/=M+1/=M-1| 多态性高时提高 | 如果从相同群体抽样降低| 系统发育研究时提高 | 通过作图观察SNP多态性和固定的变化|
