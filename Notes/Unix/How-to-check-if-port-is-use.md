---
title: 如何查看端口是否被使用
author: xuzhougeng
tags: unix, web
notebook: *NIX基础
---
# 如何查看端口是否被使用

## 使用lsof

lsof(list open file)，因为Linux里面一切都是文件，因此端口也是一种特殊的文件，如果这个端口被打开了，也就是会出现一个特殊文件表示端口。

```bash
lsof -i -P -n
# -i 选择互联网地址相关
# -P: 不对端口号进行转换
# -n: 不对IP地址进行解析成域名
```

## 使用netstat

netstat显示网络连接，地址表，交互统计数据等信息，是一个和网络相关的指令

```bash
netstat -tunlp
# -t: 显示tcp
# -u: 显示 udp
# -p: 和套字节相关的程序
# -l: 仅显示监听状态的套字节
# -n/--numeric: 显示IP地址
```