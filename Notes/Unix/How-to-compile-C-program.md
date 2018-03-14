---
title: 一文解决C程序的编译问题
author: xuzhougeng
tags: unix, 生物信息
notebook: *NIX基础
---
# 一文解决C程序的编译问题

对于C源码编译，大部分人都停留在`./configure --prefix=XXX && make && make install`这一步，大部分的程序都能顺利走完这一步，然后被安装到指定的文件下，小部分的程序会因为xxx不全而出错，然后你把这个问题放到搜索引擎上，就会找到一篇博客说用`sudo apt-get/yum install xxx` 后可以解决问题，然后问题解决了。因此从某种意义上来说编译不应该构成问题，直到你没有了root权限，直到了你开始管理服务器，不能随心所欲的用管理员权限安装程序，问题才出现了。

第一个问题是服务器的**GCC版本过低**，但不能随便升级，只能另起炉灶通过旧版本编译出新的GCC安装到类似于`/opt`,`/local/usr`或者是家目录下。我一般就用如下代码

```bash
cd ~/src
wget https://mirrors.tuna.tsinghua.edu.cn/gnu/gcc/gcc-7.2.0/gcc-7.2.0.tar.gz
tar xf gcc-7.2.0.tar.gz
cd gcc-7.2.0
./contrib/download_prerequisites
mkdir build && cd build
../configure --prefix=$HOME/opt/sysoft/gcc-7.2.0 --disable-multilib --enable-threads=posix
make -j 8 && make install
```

将编译好的gcc添加到环境变量中基本上就行了，基本上都没有遇到问题。所以如果你出错了，欢迎联系我，让我见识一下错误。

解决这个问题之后，让我们深入了解一下`./configure && make && make install`到底是做了那些事情。

众所周知，Linux的发行版是非常多的，那么代码开发的机器和实际运行的机器有很大概率是不同的，为了保证编译能通过，必须先要检查实际运行的机器是否符合要求。`configure`的首要目的就是检查依赖的**头文件**和**函数库**在目标机器上是否存在，如果存在则万事大吉，可以愉快的`make && make install`，如果不存在就得用到`configure`第二个功能，**调整编译行为**。当你用`./configure --help`时，就会出现一些输出，大致分为如下几类：

- 修改安装路径，如`--prefix`
- 调整安装目录，在`--prefix=PREFIX`的前提下，可以调整其中的lib，include, share等文件的位置
- 指定依赖工具的安装地址
- 影响编译行为的环境变量

其中最后一类的环境变量非常重要，因为它可以通过改变gcc参数影响到编译的不同阶段。这个时候，就得需要了解编译到底分为几个阶段了, 需要几个简单的例子说明.

我们都写过也编译过最简单的C程序，功能就是从屏幕上输出"hello world!", 将如下代码保存为"hello.c"

```c
#include <stdio.h>

int main(){
    printf("hello world!\n");
    return 0;
}
```

之后使用`gcc -o hello hello.c`就生成了可以执行的`hello`,虽然简单，但其实gcc做的任务不少. 它首先得在 **预处理,编译,汇编** 这一步找头文件"stdio.h", 默认查找路径是`/usr/include`, `/usr/local/include`, 如果这两个地方都没有找到，并且没有提供额外路径，那么它就会报错，比如我添加了一个默认路径没有的头文件，zlib.h，那么它就会报错，尽管系统里就有这一个文件。

```c
#include <stdio.h>
#include <zlib.h>
int main(){
    printf("hello world!\n");
    return 0;
}
# error
hello.c:2:19: fatal error: zlib.h: No such file or directory
 # include <zlib.h>
                   ^
compilation terminated.
```

那么如何解决呢？第一种方法就是提供绝对路径，即`</path/to/zlib.h>`, 但是既不美观，也不好迁移到其他环境; 第二种方法就是告诉gcc可以去哪些地方查找头文件，这个参数就是 `-I/path/to/where`, 那么最终编译代码就是`gcc -I/cluster/zlib-1.2.11/include/ -o hello hello.c`. 大部分`configure`都允许你使用"CPPFLAGS","CFLAGS"增加头文件所在路径，这两个环境变量和类似的"CXXFLAGS"基本是等价关系，因为Makefile的语法都是`$(CC) $(CPPFLAGS) $(CFLAGS) example.c -c -o example.o`.

对于一个大项目而言，通常是不会把代码放在一个文件里，而是有多个文件组成，比如说现在有两个c文件

```c
主函数：hello.c
# include <stdio.h>
int main(){
    printf("hello world!\n");
    reply();
    return 0;
}

次函数: reply.c
#include <stdio.h>
void reply(){
    printf("hi!\n");
}
```

那么为了编译出最终的程序，就需要分别编译然后调用ld进行链接，最终才能得到可执行的文件

```bash
gcc -c hello.c reply.c
gcc -o hello hello.o  reply.o
```

这种项目内先生成目标文件(.o)然后进行链接的编译行为其实不会对我们编译软件造成影响，对我们造成影响的是这个程序还用到了其他项目的函数库的时候，也就是需要其他项目的`xxx.so`文件。正如R语言我们用到最多的调用别人的R包，我们如果想用C计算sin(3.14/2)，应该不会想到动手写一个三角函数, 我们会去调用已有的函数。继续以之前hello.c 和 reply.c 为例，这不过别人已经编译了一个动态函数库libreply.so给我们调用。

```bash
gcc -shared -fPIC  -c reply.c
gcc -shared -fPIC -o libreply.so reply.o
```

当我们直接编译hello.c时候，不会通过链接这一步，因为我们的libreply.so并不在默认的搜索路径中

```c
$ gcc -o hello reply.c
/usr/lib/gcc/x86_64-redhat-linux/4.8.5/../../../../lib64/crt1.o: In function `_start':
(.text+0x20): undefined reference to `main'
collect2: error: ld returned 1 exit status
```

为了编译能够通过，就需要用`-L`和`-l`参数，第一个增加动态函数库搜索路径，第二个函数库的名字，也就是去掉"lib"和".so"剩下的部分。

```bash
gcc -o hello hello.c -L. -lreply
```

编译的确通过了，但是`./hello`运行的时候却会报错，而且这个错误我们可能非常的眼熟。

```bash
./hello
./hello: error while loading shared libraries: libreply.so: cannot open shared object file: No such file or directory
# 使用ldd检查动态库情况
$ ldd hello
	linux-vdso.so.1 =>  (0x00007ffc0a75f000)
	libreply.so => not found
	libc.so.6 => /lib64/libc.so.6 (0x00007f6224d36000)
	/lib64/ld-linux-x86-64.so.2 (0x00007f622510f000)
```

为什么会出现这个情况？这是因为动态函数库的使用分为两个情况，首先是编译的时候会用到，然后是具体执行的时候也会用到。因此，尽管我们通过了编译，但是由于libreply.so并不在系统的动态函数库搜索路径中（可以通过ldconfig -p 和 cat /etc/ld.so.conf查看），那么它依旧会因为找不到而出错。第一种解决方法是利用`LD_LIBRARY_PATH`增加动态库运行时路径，这个和`configure`无关。

```bash
export LD_LIBRARY_PATH="$HOME/temp:LD_LIBRARY_PATH"
./hello
```

这样做看似解决了问题，但其实会造成更大的隐患，你可以尝试一下编译glibc，然后把glibc的lib添加到`LD_LIBRARY_PATH`中，你就遇到神奇的**核心已转移**问题。这就是动态函数库冲突的结果，所以更好的而解决方法是用gcc的`Wl,-R/path/to/where`参数控制ld的链接行为，软连接到所需动态库上，既可以是绝对路径，也可以是相对路径。你可以用ldd去看看bioconda安装的程序的动态库位置。

```bash
gcc -L. -lreply -Wl,-R./ hello.c -o hello
$ ldd hello
	linux-vdso.so.1 =>  (0x00007ffe6cb11000)
	libreply.so => ./libreply.so (0x00007fd7609e8000)
	libc.so.6 => /lib64/libc.so.6 (0x00007fd76060d000)
	/lib64/ld-linux-x86-64.so.2 (0x000055fb9fc84000)
```

`./configure`和动态库相关的环境变量就是 **LDFLAGS**， 可以用`export LDFLAGS=-L/path/to/lib -Wl,-R/paht/to/lib`调整行为。

总结一下：使用`./configure && make && make install`编译程序的而时候，如果遇到出错，基本上就是头文件和动态库找不到，可以通过`CFLAGS/CXXFLAGS/CPPFLAGS`增加头文件路径，`LDFLAGS`增加动态库编译和搜索时的路径。不太推荐使用`LD_LIBRARY_PATH`, 会导致核心已转移。 如果程序安装说明里面没有`./configure`这一步，那么可以修改对应Makefile里的变量来解决出错。
