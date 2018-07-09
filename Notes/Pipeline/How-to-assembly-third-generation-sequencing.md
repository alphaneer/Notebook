# 如何利用三代测序组装基因组

组装算法主要分为两种，de Brujin graph(DBG)和Overlap-Layout-Consensus(OLC). DBG是用于二代测序的短序列组装，OLC最早用于Sanger测序得到的长片段序列，其中Celera是OLC算法中的代表性工具。

三代测序能够产生较长的片段，但是错误率高，而OLC算法要求两条read之间至少要有21~17个碱基是完全相同，所以无法直接套用OLC的工具。目前的三代测序得到的数据在拼接之前都要有一步纠错的过程，这一步使用诸如BLASR, MHAP, DALIGNER等三代数据比对工具寻找相似的序列，然后序列之间相互纠错。之后，使用OLC组装工具根据纠错的序列进行组装。最后组装完的结果使用二代或和三代数据进行最后的改进。

目前也有一些算法跳过原始数据一步，直接根据高噪声序列相比比对结果进行组装，最后分别利用二代或和三代数据进行纠错, 这类算法适用于重复程度比较低的小基因组，复杂基因组就不要考虑了。

## 三代组装实战

数据来自于发表在 Nature Communication 的一篇文章 "High contiguity Arabidopsis thaliana genome assembly with a single nanopore flow cell"。 这篇文章提供了 _Arabidopsis thaliana_ KBS-Mac-74 的30X短片段文库二代测序、PacBio和Nanopore的三代测序以及Bionano测序数据, 由于拟南芥的基因组被认为是植物中的金标准，因此文章提供的数据适合非常适合用于练习。根据文章提供的项目编号"PRJEB21270", 在European Nucleotide Archive上找到下载地址。

![ENA搜索](http://oex750gzt.bkt.clouddn.com/18-6-21/62677647.jpg)

```bash
## PacBio Sequal
wget -c -q ftp://ftp.sra.ebi.ac.uk/vol1/ERA111/ERA1116568/bam/pb.bam
## MinION
wget -c -q ftp://ftp.sra.ebi.ac.uk/vol1/ERA111/ERA1116595/fastq/ont.fq.gz
# Illuminia MiSeq
wget -c -q ftp://ftp.sra.ebi.ac.uk/vol1/ERA111/ERA1116569/fastq/il_1.fq.gz
wget -c -q ftp://ftp.sra.ebi.ac.uk/vol1/ERA111/ERA1116569/fastq/il_2.fq.gz
```

下载的PacBio数据以BAM格式存储，可以通过安装PacBio的smrtlink工具套装，使用其中的bam2fasta工具进行转换

```bash
# build index for convert
~/opt/biosoft/smrtlink/smrtcmds/bin/pbindex pb.bam &
# convert bam to fasta
~/opt/biosoft/smrtlink/smrtcmds/bin/bam2fasta -o pb pb.bam &
```

> PacBio的smrtlink工具套装大小为1.4G，不但下载速度慢，安装也要手动确认各种我不清楚的选项, 唯一好处就是工具很全。

### 使用Falcon进行组装

Falcon是PacBio公司开发的用于自家SMRT产出数据的基因组组装工具。对于刚接触

### 使用Canu进行组装

Canu是Celera的继任者，能用于组装PacBio和Nanopore两家公司得到的测序结果。Canu流程能自动检测服务器可用资源进行分配，使用渐进K-mer权重提高运行效率，能自动估计错误率，构建稀疏图并能输出GFA(grpahical fragment assembly)格式用于可视化展示(Bandage支持该格式)。

Canu分为三个步骤，纠错，修整和组装，每一步都差不多是如下几个步骤：

- 加载read到read数据库，gkpStore
- 对k-mer进行技术，用于计算序列间的overlap
- 计算overlap
- 加载overlap到overlap数据库，OvlStore
- 根据read和overlap完成特定分析目标
  - read纠错时会从overlap中挑选一致性序列替换原始的噪声read
  - read修整时会使用overlap确定read哪些区域是高质量区域，哪些区域质量较低需要修整。最后保留单个最高质量的序列块
  - 序列组装时根据一致的overlap对序列进行编排(layout), 最后得到contig。

这三步可以分开运行，既可以用Canu纠错后结果作为其他组装软件的输入，也可以将其他软件的纠错结果作为Canu的输入，因此下面分别运行这三步,并介绍重要的参数。

首先介绍几个全局参数：genomeSize设置预估的基因组大小，这用于让Canu估计测序深度； maxThreads设置运行的最大线程数；rawErrorRate用来设置两个未纠错read之间最大期望差异碱基数；correctedErrorRate则是设置纠错后read之间最大期望差异碱基数，这个参数需要在 **组装** 时多次调整；minReadLength表示只使用大于阈值的序列，minOverlapLength表示Overlap的最小长度。提高minReadLength可以提高运行速度，增加minOverlapLength可以降低假阳性的overlap。

**第一步**：纠错。三代测序本身错误率高，使得原始数据充满了噪音。这一步就是通过序列之间的相互比较纠错得到高可信度的碱基。主要调整两个参数

- corOutCoverage: 用于控制多少数据用于纠错。比如说拟南芥是120M基因组，100X测序后得到了12G数据，如果只打算使用最长的6G数据进行纠错，那么参数就要设置为50(120m x 50)。设置一个大于测序深度的数值，例如120，表示使用所有数据。

```bash
canu -correct \
    -p ath -d pb_ath \
    Threads=10 gnuplotTested=true\
    genomeSize=120m minReadLength=2000 minOverlapLength=500\
    corOutCoverage=120 corMinCoverage=2 \
    -pacbio-raw pb.fasta.gz
```

可以将上述命令保存到shell脚本中进行运行, `nohup bash run_canu.sh 2> correct.log &`.

这一步的overlap计算用到的算法是升级版的MHAP, MinHash，步骤还是一样都是先找到最后可能重叠的两条read，然后估计重叠的质量和长度。

注: 有些服务器没有安装gnuplot, gnuplotTested=true 可以跳过检查。这一步运行时间特别久，是项目总运行时长的一半以上。

> Canu能够利用多计算节点进行计算，对于单节点主机需要设置`useGrid=false`.

**第二步**：修整。这一步可供调整的参数几乎没有

## Pilon利用二代数据对组装进行纠错

Pilon是braodinstitute开发的java软件，能利用二代数据对组装结果进行纠错。该工具使用非常简单，只需要提供FASTA和BAM文件即可。这里的FASTA是组装后的序列，BAM文件则是用二代数据比对回组装序列得到的文件。

假设目前有如下几个文件:

- assembly.fasta: 组装后的fasta文件
- read\_r1.fastq read\_r2.fastq: 二代测序结果

第一步：利用质控软件对二代数据进行过滤. 质控软件目前常用的就是trimmomatic, fastp, cutadpter. 这里的过滤会非常严格，剔除前后低于Q38的碱基，剔除总体质量低于Q38的短读，只保留100bp以上的短读。

```bash
java -Xmx16g -jar ~/opt/biosoft/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 read_r1.fastq read_r2.fastq read_clean_1.fq.gz read_trim_1.fq.gz read_clean_2.fq.gz read_trim_2.fq.gz ILLUMINACLIP:/home/wangjw/opt/biosoft/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:38 TRAILING:38 MINLEN:100 AVGQUAL:38
```

第二步：将二代数据回贴到组装的参考序列上. 目前基因组比对都是用BWA。

```bash
# 构建索引
mkdir index
bwa index -p index/assembly assembly.fasta
# 比对
bwa mem -t 8 index/assembly read_clean_r1.fastq read_clean_r2.fastq | \
    samtools sort -@ 8 > assembly.bam
samtools index assmebly.bam
```

第三步：使用Pilon对组装结果进行纠错。它能对输入的基因组做如下的提升: 单碱基差异，小的插入缺失，稍大的局部变异，填充gap，找到局部误组装(比如说额外的gap).

Pilon特别消耗内存，官方推荐是1Mb对应1G的内存，当然也不是绝对的1:1的关系, 如果一个物种有250M，你的内存就只有256G，你可以先尝试200G内存。

```bash
java -Xmx168G -jar /home/wangjw/opt/biosoft/pilon-1.22.jar  --genome assembly.fasta --frags assembly.bam --output assembly_pillon --outdir assembly_pillon --threads 36  2> pillon.log &
```

输入参数, `--genome`提供输入参考基因组。这里的`--frags` 表示输入 < 1kb 的文库BAM， 此外还可以用 `--jumps` 输入 > 1kb 的文库BAM `-unpaired` 输入非配对的BAM。

输出参数，`--output`表示输出的前缀，`--outdir`表示输出文件夹。此外还可以用`--changes` 列出发生变化的部分,以FASTA形式保存,`--vcf` 以VCF形式保存。

除了以上简单几个参数外，还有许多复杂参数控制纠错行为。比如说`--fix` 声明对参考基因组做哪方面的改进，如"snps","indels","gaps","local", 默认是"all",也就是上面提到的这四种。

其他更多信息通过`--help`进行学习。关于Pilon的输出日志见 <https://github.com/broadinstitute/pilon/wiki/Standard-Output>

> 二代纠错步骤一般只需要1次或2次, 输入数据需要经过严格质控，避免引入新的错误。

组装软件: Canu Falcon minimap/miniasm

二代Polish: rocon, pilon

基于MUMMer的输出结果找SV: <http://assemblytics.com/>

## 参考文献

- nanopore组装拟南芥: High contiguity Arabidopsis thaliana genome assembly with a single nanopore flow cell
- 不纠错组装: Minimap and miniasm: fast mapping and de novo assembly for noisy long sequences
- 三代组装软件评测: Comprehensive evaluation of non-hybrid genome assembly tools for third-generation PacBio long-read sequence data
- Bionano: <https://github.com/RyanONeil/structome>
- Canu官方教程: <http://canu.readthedocs.io/en/latest/index.html>
