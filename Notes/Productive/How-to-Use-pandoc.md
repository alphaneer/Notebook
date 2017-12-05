---
title: 如何使用Pandoc进行markdown写作
tags: markdown
notebook: 效率工具
fignos-cleveref: On
fignos-plus-name: 图
---

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

* [如何使用Pandoc进行markdown写作](#如何使用pandoc进行markdown写作)
	* [pandoc的安装](#sec)
	* [基本markdown语法](#基本markdown语法)
	* [高级语法](#sec-1)
		* [章节编号](#sec-2)
		* [图片编号和索引](#图片编号和索引)
		* [文献引用](#文献引用)

<!-- /code_chunk_output -->

# 如何使用Pandoc进行markdown写作

工作环境:

- Windows10/MacOS
- [python](http://www.python.org/)
- [pandoc](http://www.pandoc.org/)
- [Haskell Platform](https://www.haskell.org/platform/)

## pandoc的安装 {#sec:install}

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

**段落内格式**, 包括 _italic_ (`_斜体_`), **bold** (`**粗体**`), 字体~下标~(`字体~下标~`)和字体^上标^(`字体^上标^`)以及内联代码(\`code\`)。

**插入链接**: `[github](www.github.com)`,

**插入图片**:`![](path/to/image)`.

**脚注**：(`^[脚注]`)。

**列表**：列表包括有序列表和无序列表。

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
```

**数学公式**。如果单独成行，就是用两个美元(\$)开始和结束。如果是行内的数学公式，就只要一个美元符合(\$),如 $\sum_a^b$。

当然你还可以像latex那样使用数学公式。

```equation
$$\begin{array}{ccc}
x_{11} & x_{12} & x_{13}\\
x_{21} & x_{22} & x_{23}
\end{array}$$
```

至于如何写数学公式，请自行百度Latex数学公式

## 高级语法 {#sec:highLevel}

对于一般的写作而言，学完pandoc的markdown语法基本就够了，如下部分就可以直接过掉不用看了，毕竟很多人写科技论文的首选肯定是微软的word。但其实用markdown写科技文也不是不行，你只需要往前继续再走几步而已。

Markdown的诞生就是为了尽可能简化排版，因此和排版神器LaTex相比肯定是不够看的，但LaTex强大的同时也很难学，一般只有专业人士才会使用。而pandoc目前已经有许多现成的插件用来实满足科技文写作的基本需求，例如

- 对某一个章节进行引用
- 引用图片和表格
- 引用公式
- 参考文献引用

因此，最佳方案就是利用markdown先完成初稿，使用pandoc进行格式转换，进行后期调整。

### 章节编号 {#sec:numberSection}

对章节编号并引用章节号的需求一般出现在写毕业论文或者写书的时候，一般性论文说一句“如高级语法部分所说”差不多别人就知道了。不过满足这个需求也不难，只需要安装好插件，修改pandoc的参数，markdown文本中使用专门的标识符。插件安装已经在[@sec:install]完成。

先在markdown中使用专门的标识。如果没有特殊声明，pandoc会自动为每个章节或标题赋予一个标识，具体规则不需要了解，因为我们一定要显式声明，格式为`章节{#sec:标识}`。然后以`[@sec:标识]`形进行引用。注意标识只能用全字符，不能有任何的特殊字符包括中文。

之后在命令行执行如下内容，就可以自动在输出文件中进行章节编号[^1]。

[^1]: 如果不需要对某个章节编号，则为`#section{ - }`形式

```shell
pandoc --number-sections --toc --toc-depath 3 \
    --filter pandoc-crossref \
    --latex-engine=xelatex \
    -V CJKmainfont='Ping Fang SC' \
    -V mainfont='Monaco' \
    How-to-Use-pandoc.md -o output.pdf
```

上面所达到的结果就是，见[@sec:numberSection].

和章节编号引用类似，数学公式,表格和代码的索引的语法如下

````shell
# equation
$$ math $$ {#eq:label}
# table
a   b   c
--- --- ---
1   2   3
4   5   6

: Caption {#tbl:label}
# code
```{#lst:code .haskell caption="Listing caption"}
main :: IO ()
main = putStrLn "Hello World!"
```
````

引用的语法为：

```shell
# 方法1
[@eq:label1;@eq:label2;...]
[@tbl:label1;@tbl:label2;...]
[@lst:code1;@lst:code2;...]
# 方法2
@eq:label
@tbl:label
@lst:code
```

### 图片编号和索引

关于图片要专门拿出来说下，分为一图无子图，一图多子图，多图为一图这三个情况。

一图无子图，也就是一副大图，和前面一样语法为`![Caption](file.ext){#fig:label}`，引用方法为`[@fig:label1;@fig:label2;...]`

如果一张大图里面包含多个小图，和一图无子图一样的语法，只不过引用的时候就是`[@fig:labl1]A`,手动添加子图的编号。

多子图为一图就比较麻烦，因为需要在markdown里面进行拼接，会比较麻烦，要用到HTML语法：

```shell
<div id="fig:coolFig">
![caption a](coolfiga.png){#fig:cfa width=30%}
![caption b](coolfigb.png){#fig:cfb width=60%}
![caption c](coolfigb.png){#fig:cfc width=10%}

![caption d](coolfigd.png){#fig:cfd}
![caption e](coolfige.png){#fig:cfe}
![caption f](coolfigf.png){#fig:cff}

Cool figure!
</div>
```

引用方法为`[@fig:cfa]`。

### 文献引用

和文献引用相比起来，之前的引用可能都太过简单了，我看了半天[官方文档](http://pandoc.org/MANUAL.html#citations)都没有理解他想表达什么，因为没有提供案例。最后在[如何用Markdown写论文？](http://www.jianshu.com/p/b0ac7ae98100)找到了相关案例，并且可以终于说是学会了pandoc。

在markdown中进行文献管理可以简单分为三步：

第一步，利用文献管理工具生成含有引文元数据信息的文件，如BibLaTeX,BibText. 例如使用Zotero导出参考文献为BibLatex格式，并保存为referece.bib, 见[@fig:zotero]

![从zotero导出元数据](http://oex750gzt.bkt.clouddn.com/17-12-5/69160844.jpg){#fig:zotero}

bib的内容为：

```BibLaTeX
@online{center_for_history_and_new_media_zotero_????,
	title = {Zotero 快速入门指南},
	url = {http://zotero.org/support/quick_start_guide},
	author = {{Center for History and New Media}}
}
```

第二步：在Markdown中进行引用。引用的格式`[@ID]`, 其中ID是`{`紧接的内容，这里为`center_for_history_and_new_media_zotero_????`.效果为"见实例文献[@center_for_history_and_new_media_zotero_????].

第三步: pandoc渲染.主要的参数为`--filter pandoc-citeproc`进行格式转换，`--bibliography=reference.bib`指定参考文献文件，`--csl=chinese-gb7714-2005-numeric.csl`指定参考文献格式。

```shell
pandoc --number-sections --toc --toc-depath 3 \
    --filter pandoc-crossref \
    --filter pandoc-citeproc \
    --bibliography=reference.bib \
    --csl=chinese-gb7714-2005-numeric.csl\
    --latex-engine=xelatex \
    -V CJKmainfont='Ping Fang SC' \
    -V mainfont='Monaco' \
    How-to-Use-pandoc.md -o output.pdf
```