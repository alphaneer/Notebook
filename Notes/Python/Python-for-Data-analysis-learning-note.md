---
title: Python for Data Analysis 学习笔记
tags: Python
notebook: Python笔记
---

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

# Python for Data Analysis 学习笔记

## Python内置数据结构

### 列表和相应的列表函数

Python有一些内置的列表操作函数非常好用，如：

- enumerate： 遍历同时还能进行追踪(trace)
- sorted： sorted函数能够返回排序后的结果，如果是字符串，则拆成一个个字符排序
- zip：对多个列表进行配对，返回可迭代对象

案例1： 如何对zip的序列进行解包，并记录。

```Python
seq1 = ['a','b','c']
seq2 = ['one','two','thre']
for i,(a,b) in enumerate(zip(seq1,seq2)):
    print('{0}:{1},{2}'.format(i,a,b))
```

案例2: 快速unzip已经zip的列表

```Python
seq = list(zip(seq1,seq2))
s1, s2 = zip(*seq)
```

- reversed: 序列取反 `reversed(list(range(10)))`, 返回可迭代对象

### 字典

字典几乎可以说Python中最重要的一类数据结构，

## 数据读写

数据读写分为两种方式，一种是利用Python自带的方法,一种是利用numpy和pandas提供的函数。这里主要介绍的pandas提供的一系列`read_xxx`工具，用于处理不同格式的文本。使用这些函数主要注意如下4点：

- 索引：是否需要将读取的文本中的某几列当作索引列，也就是行名，以及是否需要列名
- 类型推理和格式转换：显式或者隐式设置每一列的数据格式
- 遍历：是否对大文件遍历读取
- 数据注释： 跳过一些不需要的注视

我们这里读取拟南芥的GFF文件

## 字符操作

### Python内置函数

大部分Python内置的字符对象操作函数能够处理大部分的任务，比如说利用`split`分割字符串, 然后用`strip`删去不必要的空格和换行符。

```Python
val = 'a,b, guido'
val.split(',')
pieces = [x.strip() for x in val.split(',')]
```

一个列表里的所有子字符串可以通过`join`方法合并成一个字符串:

```Python
first, second, third = pieces
'::'.join(pieces)
```

利用Python的关键词`in`可以方便的检测字符串是否包含你需要的内容。

### pandas中向量化的字符串函数

对数据清洗的时候通常意味着要对大量的字符进行规整和标准化。这类原始数据通常可能还有一些缺失值，为了处理缺失值，pandas的Series使用了`str`属性。例如，检查某一列是否都含有某个字符串。

```Python
data = {'Dave': 'dave@google.com','Steve':'steve@gamil.com','Rob':'rob@gamil','Wes':np.nan}
data = pd.Series(data)
data.isnull()
data.str.contains('gmail')
```

或者使用更为强大的正则表达式