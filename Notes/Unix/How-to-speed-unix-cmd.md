# 如何对Linux的基本命令进行并行运算

Linux的grep默认不支持多线程，所以你即便有96个核心，也只能无奈的用一个核心，当然好在有一个神器叫做parallel, 它允许以另一种方式实现map-readuce操作。

安装

```bash
wget -4 https://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2
tar xf parallel-latest.tar.bz2
cd parallel-20180422 && ./configure --prefix=$HOME/opt
make && make install
```

用法举例：
