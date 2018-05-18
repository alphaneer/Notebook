# Python文件路径和系统操作相关库

如果想用Python替代shell脚本，我认为至少要实现如下几个方面：

- 命令行参数处理
- 文件输入和输出
- 文件系统操作
- 系统管理

涉及到两个主要的标准库: os 和 sys, 前者是操作系统交互的函数集合，后者是Python解释器运行时与外界环境交互的相关函数

**命令行参数处理**分为两个方面，参数传入和参数解析。

对于一个简单的命令行工具，先用`sys.argv`读取传入参数, 手动解析即可

```python
sys.argv[0] # 参数名
filename = sys.argv[1] # 第一个参数
sys.argv[1:] # 所有其他参数
```

如果要做一个比较复杂的命令行工具，则需要用到`argparse`, 有一定学习成本

如果你的输入全部都是文件，那么Python提供了专门的标准库`fileinput`, 比如说在Python中实现`cat`

```bash
import fileinput
import sys
# iterate all the file in the sys.argv[1:]
for line in fileinput.input():
    # output the line in standard stream, likes print
    sys.stdout.write(line)
```

**文件输入和输出**：主要是构建一个文件对象(file object), 分为文本和二进制两种形式, 随后对其操作。文件对象的构建可以有以下几个途径

- 标准输入、标准输出对象、标准错误输出对象：sys.stdin, sys.stdout, sys.stderr
- open(), 打开一个文件以供输入和输出
- io.open(), 和open()等价
- os.open()创建文件描述对象, os.fdopen()通过文件描述对象创建文件对象

```python
import os
# assume there is a test.txt in current directory
fd = os.open("test.txt", os.O_RDONLY)
# create a file object
fo = os.fdopen(fd)
```

文件对象的处理方法来自于IO库

**文件系统操作**: 会用到os(底层)和shutil(封装)

|      操作描述     |   shell |   python |
|    -----------   |  ----  |  ----    |
| 获取运行脚本所在路径| pwd | os.getcwd(), os.getcwdb()|
| 列出当前路径下文件  | ls  | os.listdir(), os.scandir()
| 切换目录           |  cd  | os.chdir(".") |
| 新建文件夹         | mkdir (-p) | os.mkdir(), os.makedirs() |
| 删除空文件夹       | rmdir (-p) | os.rmdir(), os.removedirs() |
| 删除整个文件夹      | rm -rf | shutil.rmtree() |
| 删除文件           | rm     | os.remove() |
| 移动文件,重命名     | mv | os.rename(), os.renames(), os.replace|