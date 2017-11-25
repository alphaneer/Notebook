---
title: 无root解决编译时的依赖问题
tags: unix, 服务器
notebook: *NIX基础
---

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

# 无root权限下解决编译时的依赖问题

如果你拥有最高权限，如果你只管理一台服务器，那么系统自带的包管理工具就帮你解决了所有问题。但是真实世界没有那么美好，所以我花了那么久时间去学习如何从源码开始编译一个软件。

**环境**为CentOS Linux release 7.4.1708 (Core), Linux内核version 3.10.0-693.el7.x86\_64， GCC版本为4.8.5 20150623 (Red Hat 4.8.5-16) (GCC)，

## GCC安装

首先让我们利用系统原来老旧的GCC编译器编译出最新版本的gcc吧，毕竟安装软件的时候，GCC的版本一定要过最低要求。

**第一步**： 下载gcc源码

```shell
mkdir -p ~/src && cd ~/src
wget https://mirrors.tuna.tsinghua.edu.cn/gnu/gcc/gcc-7.2.0/gcc-7.2.0.tar.gz
tar -zxvf gcc-7.2.0.tar.gz && cd gcc-7.2.0
ls
```

![](http://oex750gzt.bkt.clouddn.com/17-11-25/81177905.jpg)

**第二步**， 检查系统是否已经具备前置软件, 主要是GMP，MPFR, MPC。这些软件可以到<ftp://gcc.gnu.org/pub/gcc/infrastructure/>找到，然后下载后解压缩，并移动到gcc源码文件夹下。 可以在配置的时候用`--with-gmp, --with-mpfr --with-mpc`指定具体所在路径。

```shell
cd src
# GNU Multiple precision Library
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/gmp-6.1.0.tar.bz2 \
&& tar -jxvf gmp-6.1.0.tar.bz2 && mv gmp-6.1.0 gcc-7.2.0/gmp
# isl library
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/isl-0.18.tar.bz2 \
&& tar -jxvf isl-0.18.tar.bz2 && mv isl-0.18 gcc-7.2.0/isl
# MPFR Library
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/mpfr-3.1.4.tar.bz2 \
&& tar -jxvf mpfr-3.1.4.tar.bz2 && mv mpfr-3.1.4 gcc-7.2.0/mpfr
# MPC Library
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/mpc-1.0.3.tar.gz \
&& tar -zxvf mpc-1.0.3.tar.gz && mv mpc-1.0.3 gcc-7.2.0/mpc
```

不过更加人性化的方法是在GCC源码根目录下运行`./contrib/download_prerequisites`，可以自动搞定。

**第三步**：使用`./configure`进行配置。官方**强烈**建议, 配置所在文件夹一定要和源码所在文件夹区分开，此外configure还可以配置很多参数，我的代码如下：

```shell
mkdir build && cd build
../configure\
	--prefix=$HOME/usr \ # 指定安装路径
	--disable-multilib \ # 取消32位库编译
	--enable-threads=posix \ # 使用POSIX/Unix98作为线程支持库
```

基本上这一步不会出现太多的报错，都能够顺利生成Makefile.

**第四步**： 编译. 这步有一个小技巧就是利用多核处理器进行加速，例如`make -j2` 就是双核。

这一部分很慢很慢，因为默认设置下是3个阶段的引导(3-stage bootstrap), 以保证能够编译出完整的GCC系统并且还不会出错，你可以在配置的时候用`--disable-bootstrap`进行关闭。

**第五步**： 安装。如果你编译都成功了，那么安装也不会存在问题了， `make install`.

那么我们编译的GCC和系统自带的有什么**区别**吗？

```shell
# 从头编译
$ $HOME/usr/bin/gcc -v
Using built-in specs.
COLLECT_GCC=/home/zgxu/usr/bin/gcc
COLLECT_LTO_WRAPPER=/home/zgxu/usr/libexec/gcc/x86_64-pc-linux-gnu/7.2.0/lto-wrapper
Target: x86_64-pc-linux-gnu
Configured with: ../configure --prefix=/home/zgxu/usr --disable-multilib --enable-threads=posix
Thread model: posix
gcc version 7.2.0 (GCC)
# 系统自带
$ gcc -v
Using built-in specs.
COLLECT_GCC=gcc
COLLECT_LTO_WRAPPER=/usr/libexec/gcc/x86_64-redhat-linux/4.8.5/lto-wrapper
Target: x86_64-redhat-linux
Configured with: ../configure --prefix=/usr --mandir=/usr/share/man --infodir=/usr/share/info --with-bugurl=http://bugzilla.redhat.com/bugzilla --enable-bootstrap --enable-shared --enable-threads=posix --enable-checking=release --with-system-zlib --enable-__cxa_atexit --disable-libunwind-exceptions --enable-gnu-unique-object --enable-linker-build-id --with-linker-hash-style=gnu --enable-languages=c,c++,objc,obj-c++,java,fortran,ada,go,lto --enable-plugin --enable-initfini-array --disable-libgcj --with-isl=/builddir/build/BUILD/gcc-4.8.5-20150702/obj-x86_64-redhat-linux/isl-install --with-cloog=/builddir/build/BUILD/gcc-4.8.5-20150702/obj-x86_64-redhat-linux/cloog-install --enable-gnu-indirect-function --with-tune=generic --with-arch_32=x86-64 --build=x86_64-redhat-linux
Thread model: posix
gcc version 4.8.5 20150623 (Red Hat 4.8.5-16) (GCC)
```

不谈安装路径和版本，基本上**差别**就是在配置这一步，而这些参数就需要仔细研究了。

一个**错误**: 'Link tests are not allowed after GCC\_NO\_EXECUTABLES.' 后来发现是第三步没有在独立的文件下构建Makefile.

参考资料：

- installing GCC: <https://gcc.gnu.org/install/>
- linux下编译gcc6.2.0: <https://www.cnblogs.com/oloroso/p/5984985.html>

## Linux的编译体系

编译的常规三部曲是`./configure --prefix=$HOME/usr && make && make install`，其中最重要的一步就是`configure`，它所做的任务如下

- 检查GCC版本以及是否安装了编译所需工具
- 如果需要头文件，则默认去`/usr/include`查找
- 如果涉及到动态编译库，则默认去`/usr/lib`和`/usr/lib64`查找. 注：`lib`的函数库仅用于开机时用,提供给/bin和/sbin.

那为何需要配置？配置主要解决软件开发和软件实际安装时平台不同所导致的问题，由于平台不同，开发写的C代码需要。

## 安装zsh

```shell
wget -4 -O zsh.targz https://sourceforge.net/projects/zsh/files/latest/download
```

## 安装Tmux

准备依赖库

```shell
# libevent
cd src
wget https://github.com/libevent/libevent/releases/download/release-2.1.8-stable/libevent-2.1.8-stable.tar.gz
tar -zxvf libevent-2.1.8-stable.tar.gz && cd  libevent-2.1.8
./configure prefix=$HOME/usr && make && make install
# ncurses
cd src
wget ftp://ftp.invisible-island.net/ncurses/ncurses.tar.gz
tar -zxvf ncurses.tar.gz && configure --prefix=$HOME/usr && make && make install
```

配置环境变量，这些环境变量在configure配置时，成为GCC/CC后面的参数

```shell
# ncurses的头文件居然在include/ncurses,
CPPFLAGS="-I$HOME/usr/include/ncurses"
CFLAGS="-I$HOME/usr/include"
LDFLAGS="-L$HOME/usr/lib"
```

编译tmux

```shell
mkdir -p src && cd src
git clone https://github.com/tmux/tmux.git
cd tmux
sh autogen.sh
./configure --prefix=$HOME/usr
```

用ldd可以对程序的动态函数库进行解析

```shell
[wangjw@mgt bin]$ ldd tmux
	linux-vdso.so.1 =>  (0x00007ffd680e8000)
	libutil.so.1 => /lib64/libutil.so.1 (0x00007f3832d7d000)
	libevent-2.1.so.6 => /home6/wangjw/usr/lib/libevent-2.1.so.6 (0x00007f3832b28000)
	libresolv.so.2 => /lib64/libresolv.so.2 (0x00007f383290e000)
	libc.so.6 => /lib64/libc.so.6 (0x00007f383254b000)
	libpthread.so.0 => /lib64/libpthread.so.0 (0x00007f383232e000)
	/lib64/ld-linux-x86-64.so.2 (0x00005648487fd000)
```

事实证明"LD\_LIBRARY\_PATH"并没有多大用处。

```shell
$ wget https://cran.r-project.org/src/base/R-3/R-3.4.2.tar.gz
$ tar -zxvf R-3.4.2.tar.gz  && cd R-3.4.2/
$ ls
COPYING    INSTALL      Makefile.fw  README        VERSION       config.site  configure.ac  etc  po     src    tools
ChangeLog  Makeconf.in  Makefile.in  SVN-REVISION  VERSION-NICK  configure    doc           m4   share  tests
```