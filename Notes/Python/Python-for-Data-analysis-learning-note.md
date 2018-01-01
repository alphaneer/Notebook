---
title: Python for Data Analysis 学习笔记
tags: Python
notebook: Python笔记
---

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

# Python for Data Analysis 学习笔记

## 前言与案例

### 前言

使用Python的几个原因：

- Python目前有比较强大的数据分析库： pandas, scikit-learn
- Python和C与FORTRAN的交互比较容易，可以使用C和FROTRAN的库处理线性代数，优化问题，积分和快速傅里叶变换
- Python既能用于探索性研究，设计原型机，并且也能构建生产系统。如果使用R/SAS，后续还得用Java, C#, C++实现生产系统的代码

不用Python的理由也有：

- 运行速度太慢。不过和人力相比起来，计算机成本不算啥
- Python有一个被人嫌弃的全局解释锁(global interpreter lock, GIL)，不能使用多线程，怎么玩集群！

不过目前有PyPy等项目可以加速，初学者目前不需要担心这个运行速度问题，更大的问题是你想不出数据分析方案。

必要的Python库：Numpy(Python数值计算的基石),pandas(将R的data.frame在Python中进行实现), matplotlib(Python最棒的图形库), IPython和Jupyter(效率神器), SciPy(科学计算常用库), scikit-learn(机器学习少不了它), statsmodels(Python的统计分析库)

Python最佳发行版: Anaconda

Python版本：Python3是未来趋势

Python社区和会议：pydata,, pystatsmodels

约定俗成的几行代码:

```Python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import statsmodels as sm
```

一些术语:

- **Mugne/munging/wrangling**: 通常翻译成数据清理，指的是非结构化数据和混乱的数据经过一系列的处理成为结构和整齐的数据形式
- **pseudocode**： 伪代码，就是不能直接运行的代码，主要表现思想
- **Syntactic sugar**: 语法糖， 不是新的语言特性，仅仅是将一些语法变得更加好写，好用而已。虽然大家经常会吐槽XXX不就是语法糖嘛，但是用的时候一点都不含糊。

### 案例

在正式开始使用Python进行数据分析的时候，先了解Python能做什么能够提高我们学习的信心。作者举了几个例子，用来体现数据分析的常规流程：读取、写入数据，准备数据，格式变换，建模和计算，可视化展示。

所用数据可以在<https://github.com/wesm/pydata-book>进行下载

```bash
git clone https://github.com/wesm/pydata-book.git
```

**案例**：usa.gov的短域名数据， 展示前10个时区

```Python
%matplotlib inline
# 读取数据
import json
records = [json.loads(line) for line in open(file="../pydata-book/datasets/bitly_usagov/example.txt")]
# 数据准备，清理
from pandas import DataFrame, Series
df = DataFrame(records)
clean_tz = df['tz'].fillna('Missing')
clean_tz[clean_tz == ''] = 'Unknown'
tz_counts = clean_tz.value_counts()
# 可视化
tz_counts[:10].plot(kind='barh', rot=0)
```

![](http://oex750gzt.bkt.clouddn.com/18-1-1/51728371.jpg)

## IPython常用操作

后续的主要数据分析都是利用Ipython，IPython是加强版的Python交互终端，效率神器。一般都是通过命令行`jupyter notebook`启动浏览器版。在Jupyter下，混写markdown语法和Python代码保存为pynb格式，方便让别人理解你的分析过程。

IPython用起来和普通的Python shell差不多，但是记住如下特性能够大大提高效率

![](http://oex750gzt.bkt.clouddn.com/18-1-1/27041535.jpg)

## Python基础

### Python数据类型和数据结构

#### 列表

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

#### 字典

字典几乎可以说Python中最重要的一类数据结构，

#### 集合

### 函数

科里化（currying）: 一种计算机科学术语，指的是通过已有函数，通过修改参数构造新的函数。

```Python
from functool import partial

def add_number(x,y):
    return x+y

add_five = partial(add_number,5)
# 或
add_five = lambda y: add_number(5, y)
```

## 数据读写

数据读写分为两种方式，一种是利用Python自带的方法,一种是利用numpy和pandas提供的函数。这里主要介绍的pandas提供的一系列`read_xxx`工具，用于处理不同格式的文本。使用这些函数主要注意如下4点：

- 索引：是否需要将读取的文本中的某几列当作索引列，也就是行名，以及是否需要列名
- 类型推理和格式转换：显式或者隐式设置每一列的数据格式
- 遍历：是否对大文件遍历读取
- 数据注释： 跳过一些不需要的注视

## 数据筛选



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