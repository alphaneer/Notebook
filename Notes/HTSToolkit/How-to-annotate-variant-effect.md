---
title: 对找到的变异进行注释
tags: variant
notebook: 工具笔记
---
# 对找到的变异进行注释

变异注释软件有： VEP, snpEff, ANNOVAR, VAAST2。人类就用ANNOVAR， 软件做的还比较好。但是像我这种做植物的，就用snpEff吧。这个软件支持38,000个基因组，能整合到Galaxy,GATK等，不能再良心。

## snpEff

### 安装和配置

SnpEff是java软件，所以只需要在<http://snpeff.sourceforge.net/download.html#download>下载jar包，解压缩即可。

> 同时这也意味着你需要先安装Java

```shell
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
cp snpEff_latest_core.zip ~/biosoft
cd ~/biosoft
unzip snpEff_latest_core.zip
```

解压缩软件之后，你可能需要修改`snpEff.config`文件配置。但是重点下载安装注释数据库

#### 自动配置数据库

`snpEff`和自动配置相关的子命令是`databases`和`download`。前者配合grep检索目前支持的基因组，后者根据检索到的基因组命名进行下载。

```shell
java -jar snpEff.jar databases | grep -i thaliana
java -jar snpEff.jar download Arabidopsis_thaliana
```

实际上，某些植物你可能下载不了这些注释信息，这可能是因为网络原因，也可能是你自己拼接的基因组，所以需要自己进行配置数据库。

#### 下载已经配置好的数据库

如果你仔细观察搜索数据库时的提示信息，你会发现一串神奇的地址，这个地址里面就是已经配置好的数据库的下载地址，因此我们只要自己下载并解压就ok了。

![](http://oex750gzt.bkt.clouddn.com/17-12-28/69477217.jpg)

```bash
wget http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_ENSEMBL_BFMPP_32_24.zip
unzip snpEff_v4_3_ENSEMBL_BFMPP_32_24.zip
```

挤压之后项目文件长下面这个样子

```bash
tree -d snpEff_v4_3_ENSEMBL_BFMPP_32_24
#home/
#`-- pcingola
#    `-- snpEff
#        `-- data
#            |-- Anaplasma_phagocytophilum_str_apwi1
#            |-- Anaplasma_phagocytophilum_str_cr1007
#            |-- ...
```

也就只把snpEff这个文件夹覆盖到安装的snpEff文件夹下即可。

```bash
mv home/pcingola/snpEff/ snpEff的上一级(覆盖安装)
```

如果你不相信别人的版本，那么最终方案就是自己配置了。

#### 手动配置数据库

**step1**: 准备gene.gff(基因注释), protein.fa（蛋白序列）, cds.fa（CDS）

以拟南芥为例， 可以在ENSEMBL， TAIR, Araport上下载。这三者的基因组是一样，在注释上存在比较大的差异。

```shell
# 目前在snfEff的根目录下
mkdir ./data/TAIR10
mkdir ./data/genomes
cd ./data/TAIR10
# 下载基因注释（其实还有转座子注释）
## 建议用gffread转换成gtf格式
curl ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff -o genes.gff
# 下载CDS序列
curl ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_blastsets/TAIR10_cds_20101214_updated -o cds.fa
# 下载蛋白序列
curl ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_blastsets/TAIR10_pep_20101214_updated > protein.fa
# 下载基因组序列
cd ../genomes
curl ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas -o TAIR10.fa
```

> 注：genomes下的fa文件名要和存放的注释、CDS和蛋白序列的文件夹名相同。

一件比较尴尬的事情是，直接下载基因组序列和注释文件的第一列有可能对不上，因为染色体命名存在ensemble和NCBI两套体系。

![](http://oex750gzt.bkt.clouddn.com/17-12-28/37143155.jpg)

这就需要你用`sed`或者`awk`修改了。

**step2**: 修改snfEff目录下的`snpEff.config`

在Non-standard Databases中新增一行`xxx.genome : 物种名`, 其中xxx要和之前的xxx.fa(参考序列)命名相同，也和之前创建的文件夹名相同。

![](http://oex750gzt.bkt.clouddn.com/17-12-28/8941595.jpg)

**step3**: 构建数据库

```shell
java -jar snpEff.jar build -v TAIR10 2>&1 | tee TAIR10.build
```

最后会在`data/TAIR10`下生成`snpEffectPredictor.bin`

根据运行日志，可知如下是snpEff用于读取的文件及其路径。

![](http://oex750gzt.bkt.clouddn.com/17-12-28/74088909.jpg)

[官方文档](http://snpeff.sourceforge.net/SnpEff_manual.html#databases)还提供了RefSeq和GenBank文件构建数据库的方式。

### 对VCF进行注释

注释只需要使用`ann/eff`即可，旧版eff和新版ann等价，默认使用`ann`.

```shell
# ${id}表示一个VCF文件
java -Xmx8g -jar $snpeff TAIR10 ${id} > `basename $id`;
```

结果的VCF文件中会在INFO部分增加ANN部分（eff则是EFF）, LOF.

### 注释结果

ANN一共有16个fields，输出格式和介绍如下

```bash
Annotation      : T|missense_variant|MODERATE|CCT8L2|ENSG00000198445|transcript|ENST00000359963|protein_coding|1/1|c.1406G>A|p.Gly469Glu|1666/2034|1406/1674|469/557|  |
SubField number : 1|       2        |    3   |  4   |       5       |    6     |      7        |      8       | 9 |    10   |    11     |   12    |   13    |   14  |15| 16
```

1. Allele(ALT): 和参考基因组不同的序列
1. Annotation: 按照序列本体论(SO)描述碱基变异所带来的影响
1. Putative_impact: 简单的推测该变异的影响程度(HIGH, MODERATE, LOW, MODIFIER)
1. Gene Name: 该位点所在基因，或者是"intergenic"
1. Feature type: 下一列的feature，如transcript, motif, miRNA等
1. Feature ID: 也就是编号了
1. Transcript biotype: 至少要判断变异在"Coding",还是"Noncoding"
1. Rank/total: 变异位于第几个外显子或内含子
1. HGVS.c: 利用HGVS命名系统注释DNA
1. HGVS.p: 利用HGVS命名系统注释蛋白
1. cDNA\_position / cDNA\_len: 处于cDNA的哪个碱基
1. CDS\_position / CDS\_len: 处于CDS的哪个部位
1. Protein\_position / Protein\_len: 处于氨基酸的哪个位置
1. Distance to feature: 可选，距离最近的feature的距离
1. Errors, Warnings or Information messages： 报错，警告等信息

EFF和ANN类似，见[Input & Output files](http://snpeff.sourceforge.net/SnpEff_manual.html#input)了解。