---
title: 使用conda编译软件包
tags: conda
notebook: *NIX基础
---
# 使用conda编译软件包

## 入门篇：使用conda skeleton从PyPI编译已有的Python包

首先要保证已经安装了`Miniconda`或`Anaconda`其中一个，且保证**conda**包管理工具能够被调用。

```shell
conda install conda-build
# optional
conda upgrade conda
conda upgrade conda-build
```

如果要从PyPI编译已有的Python包，那么只要如下2个指令即可：

```shell
conda skeleton pypi pyinstrument
conda build pyinstrument
```

编译的软件包存放在`path/to/anaconda3/conda-bld/linux-64/`中，命名为`pyinstrument-0.13.2-py36hbd3c204_0.tar.bz2`.最后是用`conda install --use-local pyinstrument`进行安装。

除了PyPI外，还有cpan(perl, cpan.org), cran(R, cran.r-project.org), luarocks(lurocks.org), rpm(RPM file)等多个站点能让conda提取构建包所需信息。

## 进一步：从头构建conda包

在操作之前首先了解几个概念。第一，什么是conda的package。简单的说就是你能用`conda install [包名]`进行下载的包都算conda包，它们以打包压缩的形式存放在服务器。例如将之前`pyinstrument-0.13.2-py36hbd3c204_0.tar.bz2`解压之后会得到如下文件，

```shell
./bin:
pyinstrument
./info:
about.json  files  git  hash_input_files  hash_input.json  has_prefix  index.json  paths.json  recipe  test
./info/recipe:
conda_build_config.yaml  meta.yaml  meta.yaml.template
./lib:
python3.6
./lib/python3.6:
site-packages
./lib/python3.6/site-packages:
pyinstrument  pyinstrument-0.13.2-py3.6.egg-info
```

手动构建conda包的时候最重要的就是准备info文件下的安装说明文件(recipe)，Linux是meta.yaml, macOS是build.sh, windows则是bld.sh。因此手动构建conda包的步骤就是先创建一个和软件名的文件夹`mkdir pyinstrument`，然后在该文件夹下新建一个`meta.yaml`，之后按照下表的信息修改YAML文件，提供版本信息，源文件所在位置，依赖库是什么等。

| name	| pyinstrument |
| ---   |   ---------- |
|version|	“0.13.1” (or latest from <https://github.com/joerick/pyinstrument/releases>)|
|git_rev|	v0.13.1 (or latest from <https://github.com/joerick/pyinstrument/releases>) |
|git_url|	<https://github.com/joerick/pyinstrument.git> |
|imports|	pyinstrument|
|home|	<https://github.com/joerick/pyinstrument>|
|license|	BSD |
|license_file|	LICENSE |

```yaml
package:
  name: pyinstrument
  version: "0.13.1"

source:
  git_rev: v0.13.1
  git_url: https://github.com/joerick/pyinstrument.git

requirements:
  build:
    - python
    - setuptools
  run:
    - python

test:
  imports:
    - pyinstrument

about:
  home: https://github.com/joerick/pyinstrument
  license: BSD
  license_file: LICENSE
```

编译和安装和利用skeleton没有区别, 结果都会在conda-bld文件下生成对应的压缩包。

```shell
conda build pyinstrument
conda install --use-local pyinstrument
```