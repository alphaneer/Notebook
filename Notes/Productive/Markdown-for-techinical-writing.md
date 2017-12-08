---
title: 使用markdown进行科技论文写作
author: Xu Zhougeng
date: 2017-12-6
---

阅读这篇文章需要有一定的markdown基本知识，如果没有请查阅这篇文章[Pandoc’s Markdown 語法中文翻譯](http://pages.tzengyuxio.me/pandoci/)。操作基于MacOS, 但是仅用到开源软件，且全平台通用。用到的工具为：Zotero, pandoc。还有两个额外插件

- 插件1： pandoc-crossref(haskell-platfrom)
- 插件2：zotero-better-bibtex(zotero)

# 高级语法 {#sec:highLevel}

对于一般的写作而言，学完pandoc的markdown语法基本就够了，毕竟很多人写科技论文的首选肯定是微软的word。但其实用markdown写科技文也不是不行，你只需要往前继续再走几步，学点额外的高级语法。

Markdown的诞生就是为了尽可能简化排版，因此和排版神器LaTex相比肯定是不够看的，但LaTex强大的同时也很难学，一般只有专业人士才会使用[^1]。pandoc目前已经有许多现成的插件(pandoc-crossref)用来实满足科技文写作的基本需求，例如

- 对某一个章节进行引用
- 引用图片和表格
- 引用公式
- 参考文献引用

因此，最佳方案就是利用markdown先完成初稿，使用pandoc进行格式转换，进行后期调整。

[^1]: 目前还没有遇到如此专业的人士。

# 章节编号 {#sec:numberSection}

对章节编号并引用章节号的需求一般出现在写毕业论文或者写书的时候，一般性论文说一句“如高级语法部分所说”差不多别人就知道了。不过满足这个需求也不难，只需要安装好插件，修改pandoc的参数，markdown文本中使用专门的标识符。

先在markdown中使用专门的标识。如果没有特殊声明，pandoc会自动为每个章节或标题赋予一个标识，具体规则不需要了解，因为我们一定要显式声明，格式为`章节{#sec:标识}`。然后以`[@sec:标识]`形进行引用。注意标识只能用全字符，不能有任何的特殊字符包括中文。

之后在命令行执行如下内容，就可以自动在输出文件中进行章节编号[^2]。

[^2]: 如果不需要对某个章节编号，则为`#section{ - }`形式

```bash
pandoc --number-sections --toc --toc-depath 3 \
    --filter pandoc-crossref \
    --latex-engine=xelatex \
    -V CJKmainfont='Ping Fang SC' \
    -V mainfont='Monaco' \
    How-to-Use-pandoc.md -o output.pdf
```

上面所达到的结果就是，见[@sec:numberSection].

和章节编号引用类似，数学公式,表格和代码的索引的语法如下

````bash
# 数学公式
$$ x +y = z $$ {#eq:label}
# 表格
a   b   c
--- --- ---
1   2   3
4   5   6

: Caption {#tbl:label}
# 代码
```{#lst:code .haskell caption="Listing caption"}
main :: IO ()
main = putStrLn "Hello World!"
```
````

引用的语法基本都长这样`[@xx:xxx]`，具体一点就是：

```bash
# 方法1
[@eq:label1;@eq:label2;...]
[@tbl:label1;@tbl:label2;...]
[@lst:code1;@lst:code2;...]
# 方法2
@eq:label
@tbl:label
@lst:code
```

# 图片编号和索引

关于图片要专门拿出来说下，分为一图无子图，一图多子图，多图为一图这三种情况。

一图无子图，也就是一副大图，和前面一样语法为`![Caption](file.ext){#fig:label}`，引用方法为`[@fig:label1;@fig:label2;...]`

如果一张大图里面包含多个小图，和一图无子图一样的语法，只不过引用的时候就是`[@fig:labl1]A`,手动添加子图的编号。

多子图为一图就比较麻烦，因为需要在markdown里面进行拼接，会比较麻烦，要用到HTML语法：

```html
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

# 文献引用

和文献引用相比起来，之前的引用可能都太过简单了，我看了半天[官方文档](http://pandoc.org/MANUAL.html#citations) 都没有理解他想表达什么，因为没有提供案例。最后在[如何用Markdown写论文？](http://www.jianshu.com/p/b0ac7ae98100)找到了相关案例，并且可以终于说是学会了pandoc。

在markdown中进行文献管理可以简单分为三步：

第一步，利用文献管理工具生成含有引文元数据信息的文件，如BibLaTeX,BibText. 例如使用Zotero导出参考文献为BibLatex格式，并保存为referece.bib, 见[@fig:zotero]

![从zotero导出元数据](http://oex750gzt.bkt.clouddn.com/17-12-5/69160844.jpg){#fig:zotero}

bib的内容为：

```latex
@online{center_for_history_and_new_media_zotero,
	title = {Zotero 快速入门指南},
	url = {http://zotero.org/support/quick_start_guide},
	author = {{Center for History and New Media}}
}
```

第二步：在Markdown中进行引用。引用的格式`[@ID]`, 其中ID是`{`紧接的内容，这里为`center_for_history_and_new_media_zotero`.效果为"见实例文献[@center_for_history_and_new_media_zotero].

第三步: pandoc渲染.主要的参数为`--filter pandoc-citeproc`进行格式转换，`--bibliography=reference.bib`指定参考文献文件，`--csl=chinese-gb7714-2005-numeric.csl`指定参考文献格式。

```bash
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

等等，这是不是意味着每次都要索引的时候都需要先打开存放参考文献元数据信息的文件，然后搜索参考文献的名字，然后复制出参考文献的标识信息粘贴到文本里？如果真的是这样，我还不如直接用word写作，毕竟zotero自带的word插件是搜索起来很酷炫，如下图。

![](http://oex750gzt.bkt.clouddn.com/17-12-6/79786874.jpg){ width=50% }

好在zotero作为一个开源工具已经形成了一个生态系统，里面有各种各样的插件，其中有一个插件叫做[zotero-better-bibtex](https://github.com/retorquere/zotero-better-bibtex), 安装后新增一个标识列，如下图

![](http://oex750gzt.bkt.clouddn.com/17-12-6/28617605.jpg)

这个Citekey和原生不同，原生格式不可修改，但是这个Citekey可以在插件选项种自定义。

然后修改首选项导出的默认格式为Better BibTex Citation Key Quick Copy。该设置可以保证使用快捷复制(ctrl+shift+c/cmd+shift+c)时粘贴的数据是Citekey.

![](http://oex750gzt.bkt.clouddn.com/17-12-6/930016.jpg)

zotero-better-bibtex插件本身的选项也可以修改，这里设置为pandoc所需要的格式。
![](http://oex750gzt.bkt.clouddn.com/17-12-6/57754281.jpg)

之后在导出条目时选择'Bettet xxx'，保存为bib格式，就能和之前一样用pandoc读取。

当然还不够方便，还是需要用到鼠标，最好的方式是直接编写的时候使用一个快捷键呼出查询工具，进行选择复制格式。Atom和Sublime有对应的插件，vscode还没有。Mac平台上的Alfred3 的ZotQuery workflow能进行搜索，但是我没有搞定，可能需要额外开发了。