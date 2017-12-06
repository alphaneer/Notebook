---
title: 如何使用Pandoc进行markdown写作
tags: markdown
notebook: 效率工具
fignos-cleveref: On
fignos-plus-name: 图
---
# 如何使用Pandoc进行markdown写作

工作环境:

- 操作系统：基本上所有操作系统都行，因为后续所有软件都是开源且全平台，我用的是MacOS
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

## 使用pandoc转换markdown为其他格式