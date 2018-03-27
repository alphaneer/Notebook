---
title: Python3中文本和字节序列问题
tags: Python
notebook: Python笔记
---

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

# Python3中文本和字节序列问题

> 人类使用文本，计算机使用字节序列 ---Esther Nam 和 Travis Fisher

直到最近读取Fasta/q格式文本的时候，我想把一个fa`>`表示的序列以一行的方式读取的时候，才开始真正思考Python如何处理文本和字节的问题。

```shell
# 假设一个fasta文件里面有如下内容
> a
ATGTGT
GTGTAA
> b
GTGD
# 我的目标是Python读取的时候，每次读取的内容为`>a\nATGTGT\nGTGTAA\n`
```

如果一个你一生只需要和ASCII文本打交道，而不需要处理不同国家的文字，也不需要处理emoji，那么你将不会遇到如下问题

```Python
Traceback (most recent call last):
    File "<stdin>", line1, in <module>
UnicodeEncodeError: 'ascii' code can't encode character '\xe9' in position 3: ordinal not in range(128)
```

为了解决这个问题， 就需要了解一些更加底层的内容，首先从字符、码位和字节序列开始吧.

## 什么是字符

目前Pythn3的“字符”定义就是“Unicode字符”。Unicode标准明确区分了字符标识和具体字节表述。打个比方而言，每个人都有一套独一无二的DNA序列（字符标示），但是形成身体的不同部分则是DNA序列的具体表示方式。

- 字符的标示，也称之为**码位**（范围是0～1,114,111，相当于身份证号码)表示为4～6个十六进制且加上前缀’U+'.比如说'A'的码位就是'U+0041'
- 字符的具体表述，就是字符的码位**编码**成字节序列。比如说，A(U+0041)在UTF-8下就编码成单个字节'\x41'，而在UTF-16LE则是两个字节'\x41\x00'

从码位到字节序列是编码，而字节序列到码位则是解码。下面举例说明

```Python
s = 'baozou™'
len(s) # s有7个unicode字符
b = s.encode("utf-8") # 用UTF-8编码
b # b'baozou\xe2\x84\xa2' 9个字节
b.decode("utf-8") # 解码成人类可读
'baozou™'
```

机器存储和传输都是使用字节序列，因此需要编码器对其进行编码。人类阅读的时候再用对应的编码器进行解码。如果用错了编码器，就会出现乱码。想像成谍战片中，传输情报需要先编码，这样即便被敌人捕获了敌军也无法知道具体情报。

## 什么是字节

字节（byte）是一种计量单位，表示数据量大小，1 byte(字节) = 8 bit(比特)，其中一个比特是一个二进制的0或1， 因此1个字节在在十进制到范围是0～255，ascii字符也是256个，因此一个字节就能表示一个ascii字符。

Python3有两种二进制序列类型：不可变的bytes类型和可变的bytearray类型。bytes或bytearray对象的各个元素都是介于0～255之间的整数，即每个元素都是一个字节。虽然二进制序列其实是整数序列，但是它们的字面量表示法表明其中由ASCII文本，因此各个字节的值可能会有三种不同的方法的显示

- 可打印在ASCII范围内的字节（从空格到～），使用ASCII字符本身。
- 制表符、换行符、回车符和\对应的字节，使用转义序列`\t`,`\n`,`\r`,`\\`.
- 其他字节的值，使用十六进制转义序列，如`\x00`

这就解释前面的'baozou\xe2\x84\xa2'中我们只看到'baozou'，而没有看到'™'。

字节对象除了有和字符对象相同的方法，如`startswith,replace`等，还有自带的`decode`方法用于解码。

### 处理二进制数据的标准库

和二进制数据操作有关的标准库有两个

- struct: 主要功能是翻译打包的字节序列，比如说转换成不同类型字段组成的元祖
- codecs: 解码器，负责字节序列和码位之间的转换，Python自带了100多种解码器。

struct用于自定义二进制数据的解析方式，太过复杂，此处略去。主要介绍一下`codecs`如何使用编解码器以及常用的编码。

目前常用的编码体系如下：

- latin1: 其他编码的基础，例如cp1252,Unicode
- cp1252: 微软制定的latin1超集，也就是增加了一些有用符号
- cp437: IBM PC指定的字符集，和latin1不兼容
- gb2312: 用于编码简体中文，老古董了
- utf-8: web开发最常见的8位编码格式，与ASCII兼容。
- utf-16le: 16位编码格式

感受同一个字符串在一下不同编码体系下表现方式

```shell
for codec in ['latin_1', 'utf_8','utf_16']:
    print(codec,'hello'.encode(codec), sep="\t")
# latin_1	b'hello'
# utf_8	b'hello'
# utf_16	b'\xff\xfeh\x00e\x00l\x00l\x00o\x00'
```