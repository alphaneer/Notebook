# 使用Trinity组装转录组

## 安装Trinity

安装Trinity分为两个部分，一个是程序本身，一个程序运行要调用额外程序。程序需要通过简单的编译

```bash
wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.6.5.tar.gz
tar -xf Trinity-v2.6.5.tar.gz
cd trinityrnaseq-Trinity-v2.6.5/
make
make plugins
```

额外的必须工具官方说有三个，bowtie2，jellyfish和salmon，实际运行还会调用python2的numpy库。

```bash
# bowtie2, 解压即可
wget https://superb-dca2.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4.1/bowtie2-2.3.4.1-linux-x86_64.zip
unzip bowtie2-2.3.4.1-linux-x86_64.zip
# jellyfish, 需要编译
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.7/jellyfish-2.2.7.tar.gz
tar xf jellyfish-2.2.7.tar.gz
cd jellyfish-2.2.7
./configure --prefix=$HOME/biosoft/jellyfish-2.27
make -j 4 && make install
# salmon, 二进制解压缩
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.9.1/Salmon-0.9.1_linux_x86_64.tar.gz
unzip Salmon-0.9.1_linux_x86_64.tar.gz
mv Salmon-latest_linux_x86_64/ ~/biosoft/salmon-0.9.1
```

这些工具都**必须**添加到环境变量中，trinity才能调用。其次，还需要用pip安装numpy,即`pip install numpy`，不然运行会报错。

## 运行Trinity

官方提供了测试数据，可从<https://github.com/trinityrnaseq/trinityrnaseq/releases>找到"Trinity-v2.6.5.wExtSampleData.tar.gz"进行下载，大小为242MB.

```bash
tar xf Trinity-v2.6.5.wExtSampleData.tar.gz
```