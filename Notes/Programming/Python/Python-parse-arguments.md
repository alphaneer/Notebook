---
title: Python的命令行参数解析
author: xuzhougeng
tags: Python
notebook: Python笔记
---
# 如何使用Python对参数进行解析

主要用到一个库`argparse`，用`ArgumentParser`创建参数解析对象`ArgumentParser`，而用`add_argument()`在解析对象里添加要解析的参数. 最后用`parse_args()`进行解析，返回参数所在的命名空间。如果涉及到子命令，则还需要`add_subparsers`。

第一步：创建`ArgumentParser`对象。

```Python
import argparse
parser = argparse.ArgumentParser(description='covert all-sites vcf to fa')
```

尽管还有其他许多参数，诸如`prog`,`usage`,`add_abbre`等，但大部分情况只需要用到`description`参数，用来说明这个命令行工具的用途即可。

第二步：添加需要解析的参数

```Python
parser.add_argument('--filepath', '-f',nargs=1,required=True,help='provide a vcf file path')
```

这里表明，需要提供一个文件路径，且是必须。

第三步：解析命令

```Python
args = parser.parse_args()
```

综上，整合我已经写了的vcf转换成fa的代码，最后如下

```Python
import argparse
import re

# parse the args
parser = argparse.ArgumentParser(description='convert all-sites vcf to fa.')
parser.add_argument('--filepath','-f', nargs=1, required=True, help='a vcf file path')
args = parser.parse_args()

vcf = open(args.filepath[0])

pattern = re.compile('.*?DP=(\\d+);.*?')
current_pos = 0
min_depth = 15
seq_arr = [i for i in range(13124)]


for line in vcf.readlines():
    cols = line.split('\t')
    # get the current position
    current_chr = cols[0]
    pos = int(cols[1]) - 1
    # get the reference base and alternative base
    ref_base = cols[3]
    alt_base = cols[4]
    depth = int(re.findall(pattern = pattern, string= cols[7])[0])
    if depth > min_depth:
        if alt_base == '.':
            seq_arr[pos] = ref_base
        else:
            seq_arr[pos] = ''
    else:
        seq_arr[pos] = ''

fa = ''.join(seq_arr)
vcf.close()

with open('result.fa','wb') as f:
    f.write(fa)
```

代码还有继续的优化的余地。不过能用就行了。