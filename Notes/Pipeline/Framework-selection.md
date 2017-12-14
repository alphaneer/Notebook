---
title: 框架工具选择
author: Zhougeng Xu
date: December 5, 2017
tags: 框架, framework
notebook: 工具笔记
output: 
  pdf_document:
    toc: true
    toc_depth: 2
    highlight: tango
CJKmainfont: PingFang SC
mainfont: Monaco
---
# 框架工具选择

A review of bioinformatics pipeline framework 的作者对已有的工具进行很好的分类

![工具推荐](https://vip.biotrainee.com/assets/images/5-b9Y64zX54qeWQQir.png)

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

# nexflow起步

nextflow基于JAVA, 安装有两种方式：

- `curl -s https://get.nextflow.io | bash`
- `conda install -c bioconda nextflow`

安装之后还需要用`nexflow run hello`测试是否安装成功，安装成功后编写第一个流程`tutorial.nf`

```java
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
    file x from letters

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

```bash
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

```java
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

```bash
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

```bash
nexflow run tutorial.nf --str 'Hola mundo'
```

通过如上内容我们了解nextflow的基本工作方式，但是如何用nextflow编写流程则由后续讲解。

## nexflow基本概念：

**processes和channels**: nextflow工作流程由多个完成特定目标的process组成，每个process相互独立，互不干扰，仅仅通过异步FIFO，也就是channels进行数据交流。

**Execution abstraction**: nextflow脚本不限于运行环境，会根据实际的电脑配置运行，抽象了执行层。

**流程语言**: nextflow 有一套自己的语言定义，是[Groovy](https://en.wikipedia.org/wiki/Groovy_(programming_language))的超集， 熟悉JAVA, C/C++就能快速上手。

**配置选项**: 流程的配置文件位于当前目录下的`nextflow.config`，用于设置环境变量等参数，如

```java
env {
    PATH="$PWD/bowtie2:$PWD/tophat2:$PATH"
}
```

如此，bowtie2和tophat2就能在流程文件中直接调用

## 流程语言（pipeline language)

nextflow本身有一套语法格式，是groovy的超集，可在process和channel之间使用nextflow特有的语言特性。基本语法为：

- 屏幕输出:`println "hello world"`
- 变量赋值：`x = 1`，多变量赋值:`(a, b, c)=[10,20,'FOO']`， 和Python一样。
- 数据类型：数值，字符，布尔，日期
- 数据格式：列表（list），哈希表（maps）对应python的列表（list）和字典（dict）
    - 列表：`myList=[1,2,3,4]`
    - Maps: `scores = [ "Brett":100, "Pete":"Did not finish", "Andrew":86.87934 ]`
- 条件语句：if-else，和C语言相通

```java
x = Math.random()
if (x < 0.5) {
   println "You lost"
}
else {
    println "You win"
}
```

- 字符串:

```java
a = "world"
print "hello" + a + "\n"
println "hello $a \n"
```

列表中的字符串可以通过join方法连接。字符串支持多行，用三个引号

- 闭包，简单的说就是一组能被当作参数传递给函数的代码.

```java
square = { it * it }
[1,2,3,4].collect(square)
```

由大括号`{}`包围的表达式会被脚本解释成代码。默认闭包接受一个参数，并将其赋值给变量`it`，可以在创建时自定义。

```java
printMapClosure = { key, value ->
    println "$key = $value"
    }
[ "Yue" : "Wu", "Mark" : "Williams", "sudha" : "Kumari" ].each(printMapClosure)
```

正则表达式
文件读写：操作之前需要创建文件系统对象（file system object)

```java
myFile = file('some/path/to/my_file.file')
```

文件操作包括： 读取，创建目录，创建连接，拷贝文件，移动文件，重命名，删除，确认属性，

## processes

process分为五个部分

```java
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

channel是一个连接两个process的非阻塞无向FIFO队列，有两个特性：

1. 发出信息，异步操作，瞬间完成
1. 接受数据，阻塞操作

### Channels创建

channel可以隐式由process输出创建，或者通过_channel factory_显式定义

- create: 新建一个_channel_，`channelObj = Channel.create()`
- value：新建仅有一个值（空也是一个值）的_channel_: `expl2 = Channel.value( 'Hello there' )`
- from: 根据已有的值进行进行创建，`ch = Channel.from( 1, 3, 5, 7 )`
- fromPath: 根据文件路径创建_channel_, `myFileChannel = Channel.fromPath( '/data/some/bigfile.txt',glob: 'true' )` 支持glob匹配，不检查文件是否为空，允许有如下参数

参数名|参数描述
----|----
glob |   When true interprets characters *, ?, [] and {} as glob wildcards, otherwise handles them as normal characters (default: true)
type|	Type of paths returned, either file, dir or any (default: file)
hidden|	When true includes hidden files in the resulting paths (default: false)
maxDepth|	Maximum number of directory levels to visit (default: no limit)
followLinks|	When true it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: true)
relative|	When true returned paths are relative to the top-most common directory (default: false)

- fromFilePairs: 创建配对reads的_channel_ `Channel.fromFilePairs('/my/data/SRR*_{1,2}.fastq')`

参数名|参数描述
----| ----
type|	Type of paths returned, either file, dir or any (default: file)
hidden|	When true includes hidden files in the resulting paths (default: false)
maxDepth|	Maximum number of directory levels to visit (default: no limit)
followLinks|	When true it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: true)
size|	Defines the number of files each emitted item is expected to hold (default: 2). Set to -1 for any.
flat|	When true the matching files are produced as sole elements in the emitted tuples (default: false).

- watchPath: 时刻检查文件路径是否出现特定文件或者特定文件发生改变。

```java
Channel
   .watchPath( '/path/*.fa' )
   .subscribe { println "Fasta file: $it" }
```

当路径下出现fa的时候输出"Fasta File: xxx.fa"信息。

### Channels赋值

Channel可以通过两种方式进行复制，一种是`bind()`方法，另一种则是`<<`操作符。

```java
myChannle = Channel.create()
# 方法1
myChannel.bind( 'hello world' )
# 方法2
myChannel << 'hello world!'
```

### Channel事件监测

`subscribe()`使得每次从原始Channel输出的值都能以自定义的函数处理后的形式输出

```java
// 定义原始channel的值
source = Channel.from ( 'alpha', 'beta', 'delta' )
// 自定义函数输出
source.subscribe {  println "Got: $it"  }
// 或者类似管道形式
Channel
    .from( 'alpha', 'beta', 'lambda' )
    .subscribe { String str ->
        println "Got: ${str}; len: ${str.size()}"
     }
```

`subscribe()`还允许多个事件处理器，处理不同情形，`onNext`,`onComplete`,`onError`.

### _Channel_中的数据处理：Operators

还能使用Operator转换channel传输中的数据：过滤，塑形，分割，结合，复制，数学操作和其他。数据操作部分在进阶会继续介绍，目前只要了解其他的几个函数：`set,`ifEmpty`,`print`,`println`,`view`,`close`.

- `set`： 将Channnel的值赋予其他Channel.
- `ifEmpty`的参数值可以是闭包，在这种情况下，如果值为空则会输出闭包内结果。
- `print`和`println`两者都会输出到console的标准输出，但是后者还会在每一个输出后进行换行。`view`默认情况下类似于`println`，当设置`newLine: false`时类似于于`print`。和前两者最大的不同在于，前两者把结果输出到标准输出后结束Channel,而view会返回一个跟之前一摸一样的Channel,因此能和其他连用。
- 最后的`close`则是负责**关闭**Channel, 一般Channel会自动被Nextflow关闭，所以不太需要特别声明。

## 实例：开发出多样本SNP/InDel检测流程

核心思想：按照流程逐步开发，不断优化。

基本框架为：数据预处理（去接头，去低质量碱基）-> 建立比对索引 -> 序列比对 -> BAM排序 -> （额外处理 ->) 标记重复 -> 第一轮HC检测变异 -> 基于第一轮结果的SNP，第二轮HC 检测变异

### 配置参考基因组和FASTQ文件

问题1：如何读取单个和多个外部的文件？

解决方法：nextflow的脚本语言种提供了`file`用于文件读取，并且允许通配符如`*,?[],{}`.关于通配符，见<https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob>

```java
params.genome = "$baseDir/../ref/Athalina.fa"
```

结果会得到包含多个文件的列表。

问题2：如何区分不同样本的PE输入？

解决方法：使用_Channel_的_fromFileParis_方法

```java
params.reads = "$baseDir/../raw_data/*_{1,2}.{fq,fastq,fq.gz,fastq.gz}"
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Unable to find any reads matching: ${params.reads}" }
    .set { read_pairs }
```

最后输出新的_channel_，read_pairs，包含不同的组。

### 对参考基因组建立索引

变异检测需要对参考基因组建立两类索引：BWA索引，GATK索引。前者用于比对，后者用在HC变异检测。

要求：判断是否已经存在索引文件，如果存在则不创建。

解决方案：使用条件语句判断是否存在索引，已存在索引时使用channel,不存在时用process进行创建。

```java
# 以GATK index 为例
## 使用正则表达式构建fa_dict的文件名
fa_dict = params.genome - ~/\.fa|fasta/ + '.dict'
## 判断gatk dict是否存在，不存在时以process创建，存在时从channel中读取。
if (!file(fa_dict).exists()){
    println "building index from GATK HaplotypeCaller from ${params.genome} "
    process buildGATKDict {
        publishDir "$baseDir/../ref", mode: 'copy', overwrite: true

        input:
        file genome from genome_file

        output:
        file "${genome.baseName}.dict" into gatk_index

        """
        gatk-launch CreateSequenceDictionary -R $genome -O ${genome.baseName}.dict
        """
    }
} else{
    gatk_index = Channel.fromPath(fa_dict).view()
}
```

从中学到的一些知识点：

- .ifEmpty的闭包无法使用process.
- 不能存在同名的channel. 企图用channel和process构建同名输出channel得到的教训
- 可以用正则表达式去掉字符串的部分内容`- ~//`
- 外部条件语句可以控制Channel

### 序列比对

要求：线程和内存和分配，

最终的代码为：

```java
params.genome = "$baseDir/../ref/Athalina.fa"
params.reads = "$baseDir/../data/raw_data/*_{1,2}.{fq,fastq,fq.gz,fastq.gz}"

if (! file(params.genome).exists){
    exit 1, "${params.genome} is not exists"
}



Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Unable to find any reads matching: ${params.reads}" }
    .set { read_pairs }
```