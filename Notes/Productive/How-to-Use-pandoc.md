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
- [python](http://www.python.org/)
- [pandoc](http://www.pandoc.org/)
- [Haskell Platform](https://www.haskell.org/platform/)

## pandoc的安装{#sec:install}

最快的方法就是装一个Anaconda3, 这样子pandoc和Python都有了。至于Haskell源码编译就留待以后折腾。

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

## 高级语法 {#sec:highLevel}

对于一般的写作而言，学完pandoc的markdown语法基本就够了，如下部分就可以直接过掉不用看了，毕竟你写科技论文的首选肯定是微软的word。但其实用markdown写科技文也不是不行，你只需要往前继续再走几步而已。Markdown的诞生就是为了尽可能简化排版，因此和排版神器LaTex相比肯定是不够看的，然而LaTex不但强大而且难学，一般只有专业人士才会使用。不过pandoc目前也有许多现成的插件用来实满足科技文写作的基本需求，例如

- 对某一个章节进行引用
- 引用图片和表格
- 引用公式
- 参考文献引用

因此，最佳方案就是利用markdown先完成初稿，使用pandoc进行格式转换，进行后期调整。

### 章节编号{#sec:numer}

对章节编号并引用章节号的需求一般出现在写毕业论文或者写书的时候，一般性论文说一句“如高级语法部分所说”差不多别人就知道了。不过满足这个需求也不难，只需要安装好插件，修改pandoc的参数，markdown文本中使用专门的标识符。插件安装已经在[pandoc]{#install}完成。

先在markdown中使用专门的标识。如果没有特殊声明，pandoc会自动为每个章节或标题赋予一个标识，具体规则不需要了解，因为我们一定要显式声明，格式为`章节{#sec:标识}`。然后以`[@sec:标识]`形进行引用。注意标识只能用全字符，不能有任何的特殊字符包括中文。

之后在命令行执行如下内容，就可以自动在输出文件中进行章节编号[^1]。

```shell
pandoc -number-sections -toc --toc-depath 3 -filter pandoc-crossref How-to-Use-pandoc.md -o output.pdf
```

上面所达到的结果就是，见[@sec:numberSection].

[^1]: 如果不需要对某个章节编号，则为`#section{ - }`形式

在markdown中进行文献管理可以简单分为三步：

第一步，利用文献管理工具生成含有引文元数据信息的文件，如BibLaTeX,BibText. 例如使用Zotero导出参考文献为BibLatex格式，并保存为referece.bib, 见{@fig:zotero}

![](http://oex750gzt.bkt.clouddn.com/17-12-5/69160844.jpg){#fig:zotero}

```shell
pandoc --filter pandoc-citeproc --bibliography=myref.bib --csl=chinese-gb7714-2005-numeric.csl demo-citation.md -o demo-citation.docx
```

<https://github.com/jgm/pandoc/wiki/Pandoc-Filters>