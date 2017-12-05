---
title: 如何使用Pandoc进行markdown写作
tags: markdown
notebook: 效率工具
fignos-cleveref: On
fignos-plus-name: 图
---

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

# 如何使用Pandoc进行markdown写作

工作环境:

- Windows10/MacOS
- [python2.7](http://www.python.org/)
- [pandoc](http://www.pandoc.org/)
- [Haskell Platform](https://www.haskell.org/platform/)

## pandoc的安装

## 基本markdown语法

在学Markdown语法之前首先得理解Markdown的设计理念，方便阅读和方便发表。

> A Markdown-formatted document should be publishable as-is, as plain text, without looking like it’s been marked up with tags or formatting instructions. – John Gruber

因此，最基本的 Markdown 语法基本用手指头都能数清，都不需要动用脚趾头。

**段落**：一般而言，文章都是多个段落。Markdown分段的标准就是各个段落之间由空白行分隔。

**多级标题**：Markdown支持6级标题（如下所示），并且pandoc的markdown语法要求标题之前空一行，`#`和标题命中间有一个空格

```markdown
# 一级标题
## 二级标题
### 三级标题
#### 四级标题
##### 五级标题 = #5 五级标题
###### 六级标题 = #6 六级标题
```

第一类是**段落内格式**, 包括 _italic_ (`_斜体_`), **bold** (`**粗体**`), 字体~下标~(`字体~下标~`)和字体^上标^(字体^上标^)以及内联代码(`code`)。

还可以插入链接`[github](www.github.com)`, 插入图片也没有问题`![](path/to/image)`. 如果你对脚注有需求(`^[脚注]`)，就是这个效果^[脚注]。

第二类是**区块格式**，比如说不同级别的标题，有序列表和无序列表，索引， 代码块用三个反引号(`)开始和结束。

```markdown
* one item
* one item
* one item
    - one item
    - one item
      + one item
      + one item

1. 第一
2. 第二
3. 第三

> 这是一个索引
```

第三类是**数学公式**。如果单独成行，就是用两个美元(\$)开始和结束。如果是行内的数学公式，就只要一个美元符合(\$),如 $\sum_a^b$。

当然你还可以像latex那样使用数学公式。

```equation
$$\begin{array}{ccc}
x_{11} & x_{12} & x_{13}\\
x_{21} & x_{22} & x_{23}
\end{array}$$
```

至于如何写数学公式，请自行百度Latex数学公式

## 高级语法 {#high-level}

对于一般的写作而言，学完pandoc的markdown语法基本就够了，如下部分就可以直接过掉不用看了，毕竟你写科技论文的首选肯定是微软的word。但其实用markdown写科技文也不是不行，你只需要往前继续再走几步而已。Markdown的诞生就是为了尽可能简化排版，因此和排版神器LaTex相比肯定是不够看的，然而LaTex不但强大而且难学，一般只有专业人士才会使用。不过pandoc目前也有许多现成的插件用来实满足科技文写作的基本需求，例如

- 对某一个章节进行引用
- 引用图片和表格
- 引用公式
- 参考文献引用

因此，最佳方案就是利用markdown先完成初稿，使用pandoc进行格式转换，进行后期调整。

### 章节引用

论文中可能会用如见[高级语法](#high-level)一节这种方式来方便跳转

在markdown中进行文献管理可以简单分为三步：

第一步，利用文献管理工具生成含有引文元数据信息的文件，如BibLaTeX,BibText. 例如使用Zotero导出参考文献为BibLatex格式，并保存为referece.bib, 见{@fig:zotero}

![](http://oex750gzt.bkt.clouddn.com/17-12-5/69160844.jpg){#fig:zotero}

```shell
pandoc --filter pandoc-citeproc --bibliography=myref.bib --csl=chinese-gb7714-2005-numeric.csl demo-citation.md -o demo-citation.docx
```

<https://github.com/jgm/pandoc/wiki/Pandoc-Filters>