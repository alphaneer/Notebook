---
title: 如何对基因组序列进行注释
tags: 基因组
notebook: 分析流程
---
# 如何对基因组序列进行注释

基因组组装完成后，或者是完成了草图，就不可避免遇到一个问题，需要对基因组序列进行注释。注释之前首先得构建基因模型，有三种策略：

- 从头注释(_de novo_ prediction)：通过已有的概率模型来预测基因结构，在预测剪切位点和UTR区准确性较低
- 同源预测(homology-based prediction)：有一些基因蛋白在相近物种间的保守型搞，所以可以使用已有的高质量近缘物种注释信息通过序列联配的方式确定外显子边界和剪切位点
- 基于转录组预测(transcriptome-based prediction)：通过物种的RNA-seq数据辅助注释，能够较为准确的确定剪切位点和外显子区域。

每一种方法都有自己的优缺点，所以最后需要用EvidenceModeler(EVM)和GLEAN工具进行整合，合并成完整的基因结构。基于可靠的基因结构，后续可才是功能注释，蛋白功能域注释，基因本体论注释，通路注释等。

那么基因注释重要吗？可以说是非常重要了，尤其是高通量测序非常便宜的现在。你可以花不到一万的价格对600M的物种进行100X的普通文库测序，然后拼接出草图。但是这个草图的价值还需要你进行注释后才能显现出来。有可能你和诺贝尔奖就差一个注释的基因组。

![基因注释的重要性](http://oex750gzt.bkt.clouddn.com/18-3-10/6944182.jpg)

## 从案例中学习套路

### 陆地棉基因组注释

文章标题为“Sequencing of allotetraploid cotton (Gossypium hirsutum L. acc. TM-1) provides a resource for fiber improvement”.

**同源注释**：从Phytozome上下载了7个植物的基因组蛋白序列(Arabidopsis thaliana, Carica papaya, Glycine max, G. raimondii, Populus trichocarpa, Theobroma cacao and Vitis vinifera), 使用 _TblastN_ 将蛋白序列比对到组装序列上，E-value的阈值为1e-5. 将不同蛋白的BLAST的hits用 _Solar_ 软件进行合并。_GeneWise_ 根据每个BLAST hit的对应基因区域预测完整的基因结构。

**从头预测**：先得构建repeat-mask genome， 在这个基础上就用 _August_, _Genescan_, _GlimmerHMM_, _Geneid_ 和 _SNAP_ 预测编码区

**转录组预测**：用Tophat将RNA-seq数据比对到组装序列上，然后用cufflinks组装转录本形成基因模型。

综上，使用 _EvidenceModeler(EVM)_ 将上面的结果组装成非冗余的基因结构。进一步根据Cscore > 0.5，peptide coverage > 0.5 和CDS overlaping with TE进行筛选。还有过滤掉超过30%编码区被Pfam或Interprot TE domain的注释的基因模型。

这些基因模型使用BLASTP进行功能注释，所用数据库为SWiss-Prot和TrEMBL.蛋白功能使用InterProScan和HMMER注释，数据库为InterPro和Pfam。GO注释则是直接雇佣InterPro和Pfam注释得到的对应entry。通路注释使用KEGG数据库。

### Cardamine hirsuta基因组注释

文章标题为“The Cardamine hirsuta genome offers insight into the evolution of morphological diversity”。

**同源注释**：使用 _Genomethreader_ 以拟南芥为剪切模型，以及PlantsGDB resource上 _Brassica rapa_ (v1.1), _A. thaliana_(TAIR10), _A. lyrata_ (v6), _tomato_ (v3.6), _poplar_ (v2) 和 _A. thaliana_ (version PUT-169), _B. napus_ (version PUT-172) EST assemblies 的完整的代表性蛋白集。

**转录本预测**： 将 _C. hirsuta_ RNA-seq数据比对到基因序列，然后用cufflinks拼接

**从头预测**：转录本预测得到的潜在蛋白编码转录本使用网页工具 _ORFpredictor_ 进行预测， 同时用 _blastx_ 和 _A. thalina_ 进行比较，选择90%序列相似度和最高5%长度差异的部分从而保证保留完整的编码框(有启动子和终止子)。 这些基因模型根据相互之间的相似度和重叠度进行聚类，高度相似(>95)从聚类中剔除，保证非冗余训练集。为了训练gene finder, 它们选随机选取了2000个位点，20%是单个外显子基因。从头预测工具为 _August_ , _GlimmerHMM_, _Geneid_ 和 _SNAP_ . 此外还用了Fgenesh+, 以双子叶特异矩阵为参数进行预测。

最后使用JIGSAW算法根据以上结果进行训练，随后再次用JIGSAW对每个基因模型计算统计学权重。

可变剪切模型则是基于苗、叶、花和果实的RNA-seq比对组装结果。

GO注释使用[AHRD流程](https://github.com/groupschoof/AHRD/)

### 小结

举的2个例子都是植物，主要是植物基因组不仅是组装，注释都是一大难题。因为植物基因组有大量的重复区，假基因，还有很多新的蛋白编码基因和非编码基因，比如说玉米基因组80%以上都是重复区域。然后当我检索这两篇文章所用工具的时候，我不经意或者说不可避免就遇到了这个网站 <http://www.plantgdb.org/> , 一个整合植物基因组学工具和资源的网站，但是这个网站似乎2年没有更新了。当然这个网站也挺不错,<http://bioservices.usd.edu/gsap.html>, 他给出了一套完整的注释流程以及每一步的输入和输出情况。

此外，2017年在《Briefings in Bioinformatics》发表的"Plant genome and transcriptome annotations: from misconceptions to simple solution" 则是从五个角度对植物基因组注释做了很完整的总结

- 植物科学的常见本体
- 功能注释的常用数据库和资源
- 已注释的植物基因组意味着什么
- 一个自动化注释流程
- 一个参考流程图，用来说明使用公用数据库注释植物基因组/转录组的常规步骤

![注释流程图](http://oex750gzt.bkt.clouddn.com/18-1-23/36791187.jpg)

以上，通过套路我们对整个基因组注释有一个大概的了解，后续就需要通过实际操作来理解细节。

## 基因注释1,2,3

当我们谈到基因注释的时候，我们通常认为注释是指“对基因功能的描述”，比如说A基因在细胞的那个部分，通过招募B来调控C，从而引起病变。但是得到如下的基因结构也是注释的一种形式，也就是看似随机的ATCG的碱基排列中找到特殊的部分，而这些特殊的区域有着不一样的功能。

![gene structure](http://oex750gzt.bkt.clouddn.com/18-3-9/15619924.jpg)

在正式启动基因组注释项目之前，需要先检查组装是否合格，比如contig N50的长度是否大于基因的平均长度，使用BUSCO/CEGMA检查基因的完整性，如果不满足要求，可能输出结果中大部分的contig中都不存在一个完整的基因结构。当组装得到的contig符合要求时，就了基因组注释环节，这一步分为三步：基因结构预测，基因注释，可视化和质控。

基因结构预测和基因注释是不同的概念。前者通过计算方式，以从头预测、同源比对的方式得到**可能**的基因结构，相当于收集线索。后者是基于前者收集的线索进行梳理，最终确定“最可能”的模型。这一点告诉我们，面对非模式植物的注释时，一定要谨慎，不要盲目使用。也就是说我们需要通过进一步通过可视化的方式，抽样或者挑选自己的目标基因进行检查，因为从某种程度上，人类模式识别能力还是最一流的。

### 重复序列屏蔽

真核生物的基因组存在大量的重复序列，植物基因组的重复序列甚至可以高达80%。尽管重复序列对维持染色体的空间结构、基因的表达调控、遗传重组等都具有重要作用，但是却会导致BLAST的结果出现大量假阳性，增加基因结构的预测的计算压力甚至影响注释正确性。基因组中的重复按照序列特征可以分为两类：串联重复(tandem repeats)和散在重复(interspersed repeats).

![人类中的重复序列划分](http://oex750gzt.bkt.clouddn.com/18-3-14/14545377.jpg)

鉴定基因组重复区域的方法有两种：一种基于文库(library)的同源(homology)方法，该文库收集了其他物种的某一种重复的一致性序列，通过相似性来鉴定重复；另一种是从头预测(_de novo_)，将序列和自己比较或者是高频K-mer来鉴定重复。

## 功能注释

- IPR, Pfam: interproscan
- GO: interproscan
- KEGG: 本地BLAST
- uniport: 本地BLAST

## 附录

基因组注释的常用软件：

- 重复区域
  - RepeatMasker：识别基因组中的可能重复
  - RepeatModeler: 识别新的重复序列
  - LTR-FINDER: <http://tlife.fudan.edu.cn/ltr_finder/>
- 从头预测
  - Augustus
  - Fgeneesh
- 同源预测
  - GeneWise
  - Trinity
- 自动化注释
  - GLEAN
  - EvidenceModeler
- BRAKER1: 使用GeneMark-ET和AUGUSTUS基于RNA-Seq注释基因结构
- 流程
  - PASA：真核生物基因的转录本可变剪切自动化注释项目，需要提供物种的EST或RNA-seq数据
  - MAKER
- 可视化
  - IGV
  - JBrowse/GBrowse

参考文献和推荐阅读：

- NCBI真核生物基因组注释流程<https://www.ncbi.nlm.nih.gov/genome/annotation_euk/process/>
- 真核基因组注释入门: "A beginner’s guide to eukaryotic genome annotation"
- 二代测序注释流程:Comparative Gene Finding: "Annotation Pipelines for Next-Generation Sequencing Projects"
- 基因组转录组注释策略: "Plant genome and transcriptome annotations: from misconceptions to simple solution"
- 重复序列综述: "Repetitive DNA and next-generation sequencing: computational challenges and solutions"
- MAKER2教程: <http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018>
- 《生物信息学》 樊龙江: 第1-5章: 基因预测与功能注释
- 《NGS生物信息分析》 陈连福： 真核生物基因组基因注释

## 环境准备

这一部分主要准备练习用数据和安装MAKER。数据下载比较容易，工具安装则可能比较麻烦。

```bash
# Cardamine hirsutat基因组数据
mkdir chi_annotation && cd chi_annotation
wget http://chi.mpipz.mpg.de/download/sequences/chi_v1.fa
cat chi_v1.fa | tr 'atcg' 'ATCG' > chi_unmasked.fa
# Cardamine hirsutat转录组数据
wget -4 -q -A '*.fastq.gz' -np -nd -r 2 http://chi.mpipz.mpg.de/download/fruit_rnaseq/cardamine_hirsuta/ &
wget -4 -q -A '*.fastq.gz' -np -nd -r 2 http://chi.mpipz.mpg.de/download/leaf_rnaseq/cardamine_hirsuta/ &
# maker2教程数据
wget http://weatherby.genetics.utah.edu/data/maker_tutorial.tgz
tar xf maker_tutorial
```

MAKER依赖的软件包比较多，最容易的方法就是用bioconda安装，你只需要等待半小时就似乎得到了一个可用的环境。下面部分仅提供另一种折腾思路，不推荐模仿。

这里安装的MAKER-P，包含标准的MAKER，以及植物注释的相关工具。

首先你需要先通过<http://www.yandell-lab.org/software/maker.html>里注册链接，表明自己是用于学术用途，而不是商业。目前稳定版是2.31.9，开发版是3.01，我使用稳定版。下载之后，解压缩，按照里面的要求逐步安装。

> 如无特殊说明，软件下载后存放在`~/src`目录下，生信软件安装在`~/opt/biosoft`,系统工具安装在`~/opt/sysoft`, 单个二进制软件存在`~/opt/bin`.

**OpenMPI**: MPI(Message Passing Library)库，可以提高运行速度，选择性安装

```bash
cd ~/src
wget https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.gz
tar xf openmpi-3.0.0.tar.gz
cd openmpi-3.0.0
./configure --prefix=$HOME/opt/sysoft/openmpi-3.0.0
make && make install
# 环境变量
export PATH=$HOME/opt/sysoft/openmpi-3.0.0/bin:$PATH
export LD_RUN_PATH=$HOME/opt/sysoft/openmpi-3.0.0/lib:${LD_RUN_PATH:+:${LD_RUN_PATH}}
export LD_LIBRARY_PATH=$HOME/opt/sysoft/openmpi-3.0.0/:${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
```

**Perl模块**, 我选择编译一个最新的Perl，因为自己编译Perl的好处就在于之后的perl模块都会安装到自己的Perl目录下，而不会对系统造成影响。

```bash
cd ~/src
wget -4 http://www.cpan.org/src/5.0/perl-5.26.1.tar.gz
tar xf perl-5.26.1.tar.gz
cd perl-5.26.1
./Configure -des -Dprefix=$HOME/opt/sysoft/perl-5.26.1
make test
make install
```

为了方便后续安装，先下载一个`cpanm`, 并增加国内镜像地址

```bash
wget https://cpan.metacpan.org/authors/id/M/MI/MIYAGAWA/App-cpanminus-1.7043.tar.gz
tar xf App-cpanminus-1.7043.tar.gz
cd App-cpanminus-1.7043
perl Makefile.PL
make test && make install
echo 'alias cpanm="~/opt/sysoft/perl-5.26.1/bin/cpanm --mirror http://mirrors.163.com/cpan --mirror-only"' >>~/.bashrc
```

重启终端后安装所有Perl模块

```bash
cpanm DBI DBD::SQLite forks forks::shared File::Which Perl::Unsafe::SIgnals Bit::Vector Inline::C IO::All IO::Prompt Bundle::BioPerl Text::Soundex
```

**BLAST**，BLAST有两个版本可供选择, WuBLAST或者NCBI-BLAST，我个人倾向于NCBI-BLAST，并且推荐使用编译后二进制版本，因为编译实在是太花时间了

```bash
cd ~/src
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+-x64-linux.tar.gz
tar xf ncbi-blast-2.7.1+-x64-linux.tar.gz
mv ncbi-blast-2.7.1+ ~/opt/biosoft
# 环境变量
export PATH=~/opt/biosoft/ncbi-blast-2.7.1+/bin:$PATH
```

**SNAP**: 基因从头预测工具，最后一个版本是2013/11/29，效果还行，在处理含有长内含子上的基因组上表现欠佳

```bash
# 安装
cd ~/src
wget http://korflab.ucdavis.edu/Software/snap-2013-11-29.tar.gz
tar xf snap-2013-11-29.tar.gz
cd snap
make
cd ..
mv snap ~/opt/biosoft
# 环境变量
export Zoe=~/opt/biosoft/snap/Zoe
export PATH=~/opt/biosoft/snap:$PATH
```

**RepeatMasker**: 用于注释基因组的重复区，需要安装RMBlast, TRF，以及在<http://www.girinst.org>注册以下载Repbase

```bash
# trf
cd ~/src
wget http://tandem.bu.edu/trf/downloadstrf409.linux64
mv trf409.linux64 ~/opt/bin/trf
chmod a+x ~/opt/bin/trf
# RMBlast
cd ~/src
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-src.tar.gz
wget http://www.repeatmasker.org/isb-2.6.0+-changes-vers2.patch.gz
tar xf ncbi-blast-2.6.0+-src
gunzip isb-2.6.0+-changes-vers2.patch.gz
cd ncbi-blast-2.6.0+-src
patch -p1 < ../isb-2.6.0+-changes-vers2.patch
cd c++
./configure --with-mt --prefix=~/opt/biosoft/rmblast --without-debug && make && make install
# RepeatMasker
cd ~/src
wget http://repeatmasker.org/RepeatMasker-open-4-0-7.tar.gz
tar xf RepeatMasker-open-4-0-7.tar.gz
mv RepeatMasker ~/opt/biosoft/
cd ~/opt/biosoft/RepeatMasker
## 解压repbase数据到Libraries下
## 配置RepatMasker
perl ./configure
## 安装成功后的输出
Congratulations!  RepeatMasker is now ready to use.
The program is installed with a the following repeat libraries:
  Dfam database version Dfam_2.0
  RepeatMasker Combined Database: Dfam_Consensus-20170127, RepBase-20170127
# 添加环境变量
export PATH=~/opt/biosoft/maker:$PATH
```

**Exonerate 2.2**: 配对序列比对工具，提供二进制版本, 功能类似于GeneWise，能够将cDNA或蛋白以gao align的方式和基因组序列联配。

```bash
cd ~/src
wget http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.2.0-x86_64.tar.gz
tar xf exonerate-2.2.0-x86_64.tar.gz
mv exonerate-2.2.0-x86_64 ~/opt/biosoft/exonerate-2.2.0
# .bashrc添加环境变量
export PATH=~/opt/biosoft/exonerate-2.2.0:$PATH
```

安装MAKER本体

```bash
cd ~/src
tar xf maker-2.31.9.tgz
mv maker ~/opt/biosoft
cd ~/opt/biosoft/maker/src
perl Buidl.PL # 根据输出信息查漏补缺，选择是否使用MPI
./Build install
# 添加环境变量
export PATH= ~/opt/biosoft/maker/bin:$PATH
```