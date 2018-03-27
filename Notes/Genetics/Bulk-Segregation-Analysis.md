# bulk segregation analysis: 混池测序分析

基于BSA分析思想的一系列方法文章

- mapping-by-sequencing: "SHOREmap: simultaneous mapping and mutation identification by deep sequencing"
- whole-genome sequencing: "Caenorhabditis elegans mutant allele identification by whole-genome sequencing."
- pool-seq: "PoPoolation2: identifying differentiation between populations using sequencing of pooled DNA samples (Pool-Seq)"
- MutMap: "Genome sequencing reveals agronomically important loci in rice using MutMap"
- QTL-seq: "QTL-seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations"
- pnome: "PpeTAC1 promotes the horizontal growth of branches in peach trees and is a member of a functionally conserved gene family found in diverse plants species."
- 杂合基因组，林木:"Exploring DNA variant segregation types in pooled genome sequencing enables effective mapping of weeping trait in Malus"

所使用的策略：

- variant scarcity or density mapping
- variant discovery mapping
- SNP index
- bulk segregant linkage mapping
- delta SNP index mapping
- SNP ratio mapping(SRM)
- mutatnt allele frequency(MAF)
- allelic distance(AD)
- homozygosity mapping
- Mutant allele frequency and density (MAFD)
- Allele frequency directional difference and density(AFDDD)

开发了一系列软件: SHOREmap, CloudMap, SNPtrack, MegaMapper, MMAPPR, EXPLoRA, GIPS，每个软件的策略都有所不同，但是利用的信息差不多就是如下三个：

- 变异位点等位基因频率(variant allele frequency)
- 变异位点密度(variant density)
- 变异位点距离(variant distance)

使用F1的群体BSA分析时，后代的变异可能情况:

![图示](http://oex750gzt.bkt.clouddn.com/18-3-16/58333063.jpg)

推荐综述:

- Using next-generation sequencing to isolate mutant genes from forward genetic screens
- Bulked sample analysis in genetics, genomics and crop improvement