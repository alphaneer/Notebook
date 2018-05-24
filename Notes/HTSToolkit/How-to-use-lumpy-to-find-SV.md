# 使用LUMPY检测结构变异

LUMPY是一款基于概率框架检测结构变异(structure variants)的软件, 它根据read-pair, split-read, read-depth和其他先验知识寻找基因组上可能的结构变异。

软件在编译的时候会先安装HTSLIB，根据Makefile, 需要预先安装好curl和zlib. 此外还推荐安装Python的Pysam和Numpy，Samtools(0.1.18+)，SAMBLASTER(0.1.19+), sambamba。

安装完成之后可以用测试数据进行测试, 数据下载地址为<http://layerlab.org/lumpy/data.tar.gz>, 要存放在`/lumpy-sv/data`下