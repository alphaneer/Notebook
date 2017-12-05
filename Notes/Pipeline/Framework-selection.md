---
title: 框架工具选择
author: Zhougeng Xu
date: December 5, 2017
tags: 框架, framework
notebook: 工具笔记
<<<<<<< HEAD
output: 
  pdf_document:
    toc: true
    toc_depth: 2
    highlight: tango
CJKmainfont: PingFang SC
mainfont: Monaco
=======
>>>>>>> 986a83c90fedb3b8abd99d782db35aba8e4d7bb5
---

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

* [框架工具选择](#框架工具选择)
	* [nexflow起步](#nexflow起步)
	* [nexflow基本概念：](#nexflow基本概念)
	* [流程语言（pipeline language)](#流程语言pipeline-language)
	* [Channels](#channels)

<!-- /code_chunk_output -->

# 框架工具选择

A review of bioinformatics pipeline framework 的作者对已有的工具进行很好的分类

![image https://vip.biotrainee.com/assets/images/5-b9Y64zX54qeWQQir.png](https://vip.biotrainee.com/assets/images/5-b9Y64zX54qeWQQir.png)

作者的看法：

1. implicit，也就是Make rule语法更适合用于整合不同执行工具
1. 基于配置的流程更加稳定，也比较适合用于集群分配任务。

最后作者建议是：

- 如果实验室既不是纯粹的生物学试验（不需要workbench这种UI界面），也不需要高性能基于类的流程设计， 不太好选， 主要原则是投入和产出比
- 如果实验室进行的是**重复性**的研究，那么就需要对数据和软件进行**版本控制**， 建议是 configuration-based pipelines
- 如果实验室做的是探索性的概念证明类工作（exploratory proofs-of-concept)，那么需要的是 DSL-based pipeline。
- 如果实验室用不到高性能计算机(HPC)，只能用云服务器，就是server-based frameworks.

目前已有的流程可以在[awesome-pipeline](https://github.com/pditommaso/awesome-pipeline) 进行查找。

就目前来看，pipeline frameworks & library 这部分的框架中 [nextflow](https://github.com/nextflow-io/nextflow) 是点赞数最多的生物学相关框架。所以我就开始学习了nextflow, 只可惜nextflow在运行时需要创建fifo，而在NTFS文件系统上无法创建，所以Ubuntu On Windows10是玩不转的。

## nexflow起步

nextflow基于JAVA, 安装有两种方式：

- `curl -s https://get.nextflow.io | bash`
- `conda install -c bioconda nextflow`

安装之后还需要用`nexflow run hello`测试是否安装成功，安装成功后编写第一个流程`tutorial.nf`

```nexflow
#!/usr/bin/nextlfow

params.str = 'Hello world!'

process splitLetters {
    output:
    file 'chunk_*' into letters mode flatten

    """
    printf '${params.str}' | split -b 6 - chunk_
    """
}

process convertToUpper{
    input:
    file x fomr letters

    output:
    stdout result

    """
    cat $x | tr '[a-z]' '[A-Z]'
    """
}

result.susribe {
    println it.trim()
}
```

在命令行执行如下命令

```shell
nextflow run tutorial
# 终端输出内容
N E X T F L O W  ~  version 0.25.1
Launching `tutorial.nf` [jolly_ride] - revision: 4a11bf2927
[warm up] executor > local
[92/fd0b10] Submitted process > splitLetters
[d3/8f6c61] Submitted process > convertToUpper (2)
[2c/39e326] Submitted process > convertToUpper (1)
WORLD!
HELLO
```

nextflow在运行时会在当前路径下生成"work"目录，用于记录运行时每一步的数据。这意味着什么？这意味着你修改了脚本其中一个部分的时候，继续运行可以基于已有的数据，而不需要重新重头开始。我们可以尝试修改其中`convertToUpper`部分。

```nexflow
process convertToUpper {

    input:
    file x from letters

    output:
    stdout result

    """
    rev $x
    """
}
```

在原来执行的方式上加上`-resume`参数，

```shell
nextflow run tutorial.nf -resume
# 终端输出内容
N E X T F L O W  ~  version 0.25.1
Launching `tutorial.nf` [berserk_torvalds] - revision: 211d9426dc
[warm up] executor > local
[92/fd0b10] Cached process > splitLetters
[bc/dd002c] Submitted process > convertToUpper (2)
[9b/ef1aa9] Submitted process > convertToUpper (1)
!dlrow
olleH
```

你会发现前后两部有一个共同的部分`[92/fd0b10] Cached process > splitLetters`, 也就是说修改后的代码时直接从中间某一步继续，而不是从头到位跑下来。

nextflow还支持外部输入参数，覆盖已有的设置`params.str = 'Hello world!'`。

```shell
nexflow run tutorial.nf --str 'Hola mundo'
```

通过如上内容我们了解nextflow的基本工作方式，但是如何用nextflow编写流程则由后续讲解。

## nexflow基本概念：

**processes和channels**: nextflow工作流程由多个完成特定目标的process组成，每个process相互独立，互不干扰，仅仅通过异步FIFO，也就是channels进行数据交流。

**Execution abstraction**: nextflow脚本不限于运行环境，会根据实际的电脑配置运行，抽象了执行层。

**流程语言**: nextflow 有一套自己的语言定义，是[Groovy](https://en.wikipedia.org/wiki/Groovy_(programming_language))的超集， 熟悉JAVA, C/C++就能快速上手。

**配置选项**: 流程的配置文件位于当前目录下的`nextflow.config`，用于设置环境变量等参数，如

```nextflow
env {
    PATH="$PWD/bowtie2:$PWD/tophat2:$PATH"
}
```

如此，bowtie2和tophat2就能在流程文件中直接调用

## 流程语言（pipeline language)

在process和channel之间能够使用nextflow特有的语言特性。

输出:`println "hello world"`
变量：`x = 1`
列表： `myList=[1,2,3,4]`
Maps: `scores = [ "Brett":100, "Pete":"Did not finish", "Andrew":86.87934 ]`， 相当于Python的字典。
多变量赋值:`(a, b, c)=[10,20,'FOO']`， 和Python一样。
条件语句：

```if-else
x = Math.random()
if (x < 0.5) {
   println "You lost"
}
else {
    println "You win"
}
```
字符串:
```python
a = "world"
print "hello" + a + "\n"
println "hello $a \n"
```
列表中的字符串可以通过join方法连接。字符串支持多行，用三个引号

闭包，简单的说就是一组能被当作参数传递给函数的代码
```python
square = { it * it }
[1,2,3,4].collect(square)
```

正则表达式
文件读写：操作之前需要创建文件系统对象（file system object)
```
myFile = file('some/path/to/my_file.file')
```
文件操作包括： 读取，创建目录，创建连接，拷贝文件，移动文件，重命名，删除，确认属性，



process分为五个部分
```nextflow
process < name > {

   [ directives ]

   input:
    < process inputs >

   output:
    < process outputs >

   when:
    < condition >

   [script|shell|exec]:
   < user script to be executed >

}
```

最重要的部分是脚本，nextflow默认以BASH脚本执行

## Channels

channel是一个连接两个process的非阻塞无向FIFO队列，由两个特性：
1. 发出信息，异步操作，瞬间完成
2. 接受数据，阻塞操作

channel可以隐式由process输出创建，或者用如下方法显式定义
- create `channelObj = Channel.create()`
- empty 
- from `ch = Channel.from( 1, 3, 5, 7 )`
- fromPath `myFileChannel = Channel.fromPath( '/data/some/bigfile.txt' )` 支持glob匹配
- fromFilePairs `Channel.fromFilePairs('/my/data/SRR*_{1,2}.fastq')`
- value `expl2 = Channel.value( 'Hello there' )`
- watchPaht

还能使用Operator转换channel传输中的数据
- 过滤
- 塑形
- 分割
- 结合
- 复制
- 数学操作
- 其他
