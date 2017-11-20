---
title: Python for Data Analysis 学习笔记
tags: Python
notebook: Python笔记
---

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

# Python for Data Analysis 学习笔记

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