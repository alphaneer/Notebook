---
title: 使用make构建分析流程
tags: make, unix
notebook: *NIX基础
---
# 使用make构建分析流程

除了shell script以外，还有一种更加强大的工具--`make`, 它最早是用于解决大型软件的编译过程中的依赖问题。它是一种[Domain specific language](https://en.wikipedia.org/wiki/Domain-specific_language), 换句话说就是make专注于各种文件依赖关系，非常符合生信分析就是一个格式变成另一个格式的风格。

## 事先准备

准备时间: 1~20min

首先你需要保证你用的机器上安装了GNU Make 和 Python2/3(numpy, matplotlib), 然后下载后续操作所需要的材料

```bash
wget http://swcarpentry.github.io/make-novice/files/make-lesson.zip
unzip make-lesson.zip
cd make-lesson
```

对于Python，建议安装Anaconda/Miniconda.

## 简介：为什么是make

## 写一个简单的Makefile

```Make
<target> : <prerequisites>
[tab]  <commands>
```

## D.R.Y: 不要重复工作

自动变量: `$@`, `$<`, `$^`, `$?`

## 数据和代码依赖

## 模式规则

`%`仅能用与目标和依赖文件，而不能用在执行命令中。在执行命令中需要用`$*`指代**匹配部分**

## 变量

使用在Makefile中使用`include config.mk`实现配置分离

## 函数

- `wildcard`: 文件统配
- `patsubst`: 模式替换

## 如何做一个make自文档

使用`sed`进行处理，同时在文档中用`##`进行标注。

```make
.PHONY : help
help: Makefile
        @sed -n 's/^##//p' $
```

## 参考资料

- [Automation and Make](http://swcarpentry.github.io/make-novice/)
- [Make 命令教程](http://www.ruanyifeng.com/blog/2015/02/make.html)