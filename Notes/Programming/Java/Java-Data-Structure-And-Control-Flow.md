---
title: Java学习笔记（二）Java基本语法
tags: Java
notebook: Java笔记
---
# Java学习笔记（二）Java基本语法

## Java的数据类型

Java是一门纯粹的面向对象编程语言，除了8个基本数据类型不是对象以外，其他的一切的都是对象。那么问题来了，这8个数据类型分别是啥？

![](http://oex750gzt.bkt.clouddn.com/17-12-13/91989074.jpg)

不难发现，Java的数据类型和C语言基本上都是一致的，两者都不包含**字符串**这种数据类型。尽管Java可以使用`String str=Hello World!"`这种方法来定义字符串，但其实和`char str[] = "Hello World!"`本质是一样的，也就是通过组合字符数据类型和数组来保存字符串。因此就对单引号和双引号进行了区分，而在Python里面没有单个字符，只有长度为1的字符串。

这几种基本数据类型可以相互转换，分为强制转换和自动转换两种类型，如下代码

```Java
//强制转换
int a = (int) (3.14 * 2.13)
// a =6
// 自动转换
float b = a
// b = 6.0
```

除了基本数据类型外，Java还有一种称之为引用类型，包括类、接口和数组类型以及特殊的null类型。所谓引用类型就是对一个对象的引用，听起来好像是C语言的指针一样，其实就是，只不过Java不再使用指针这个说法。

## Java的运算符

Java的运算符有如下几种

- 算术运算符: +, -, *, \, %
- 赋值运算符: =, +=, -= ...
- 比较运算符: >, >=, <, <=, ==
- 逻辑运算符: &&, ||, !
- 位运算符: &, |, ~, ^, <<, >>, >>>. 从来没有用过，我很尴尬，只知道效率很高。
- 三目运算符:(expression) ? 条件为真时 : 条件不为真时

通过将运算符和数据结合就构成了赋值表达式，而当使用多个类型的运算符时需要注意运算符之间的优先级关系。不过真实的编程世界不会存在多个运算符同时登场让你怀疑人生，以为自己在参加什么证书考试。如果有哪个程序员闲的蛋疼真的要这样做，直接开除好了，否则留着心累吗。

将运算符和不同数据类型进行结合就是**运算表达式**，运行结果不会直接丢掉，而是保存在变量中。

## Java的变量

编程的本质可以认为就是对内存中的数据的访问和修改。计算机在内存中开辟空间用于存放数据，那么程序如何访问这些数据并进行修改呢？尽管每个数据在内存上都有一个确定的位置，但是程序员肯定是不会通过手动输入地址的方式来获取数据，而是为存放在内存中的数据赋予一个变量名，通过访问变量名的方式来对数据做一系列的修改。

Java和C都是强类型的编程语言，强类型有两种特征：1. 所有的变量都必须先声明后才能使用；2. 指定类型的变量只能接受类型与之匹配的值。也就是说下面这段Python代码在Java编译时会出错。

```Python
c = 1
c = 'abc'
```

在Java和C这种强类型语言中，你需要这样写:

```Java
int num = 1
char str[] = "abc"
```

强类型有一定的优点，能够提高运行速度，降低一些编码错误。对于我这种没有好好深入学习强类型语言的人而言，目前还在思索如何读取不定行数字符串输出成一个。不难发现Java中的变量定义由3部分组成：基本数据类型 + 变量名（声明符） + 变量值（可选）。当然后面指针的定义会稍微复杂。

## Java的三种程序结构

1996年，计算机科学家Bohm和Jacopini证明了“无论一个算法是否简单，还是足够复杂，都可以用顺序结构，选择结构和循环结构这三种基本结构组合而成”。因此所有计算机编程语言都会具备这三个基本结构，以我浅显码代码生涯，还没有遇到一门语言不存在这三个基本组合。

**顺序结构**：顺序结构就是从上到下运行，行云流水，没有犹豫（条件语句），没有彷徨（循环语句），一直到程序运行结束为止。之前的HelloWorld祷告与就是如此。非常直白，不需要多余的解释。

**分支结构**：Java提供了两种常见的分支控制结构:if语句和switch语句。Python之父觉得switch没啥必要，所以Python里只有if-else。

**循环结构**：循环有三种，while, do while,for. 并且循环可以套循环。谈及循环结构就一定要说说`break`和`continue`两者的区别，前者是彻底不干，后者是当前循环内容不做，继续搞下一个循环。

## 数据结构：数组

基本数据类型的元素在内存中随机存放，数组是一种最常见的数据结构，用于将**相同类型**的数据存放在同一块内存区。在C/C++里面，数组一般都要和指针联系在一起。Java没有说指针这个概念，提出了引用类型。

### 数组定义

Java支持两种方式定义数组:

```Java
type[] arrayName;
type arrayName[];
```

我看的书里面推荐Java里使用数组时建议以第一种格式进行。因为第一种有着更好的可读性，一看就知道int和int[]是两种不同的数据类型，int是基本类型，而int[]就是引用类型。但是在C/C++里面，似乎只能用第二种格式。为什么不能用第一种格式呢？当我看到`int *(&arry)[10] = ptrs`，我觉得这就是一个原因吧。

还有Java定义数组**不能指定数组的长度**。因为定义数组只是定义了一个引用变量，并未指向任何有效的内存空间，也就是说没有内存空间来存储数组元素。只有当数组进行初始化后才能使用这个数组。

### 数组初始化

数组初始化意味着在内存中开辟一个空间，这个空间存放数组元素，每个数组元素都需要被赋值。这个值可以手动指定，也可以系统自动赋予，或者就是null数据类型。Java初始化数组由两种方式：静态初始化和动态初始化。

- 静态初始化：程序员显式指定每个数组元素的初始值，系统决定数组长度
- 动态初始化: 程序员显式指定数组长度，系统分配初始值

首先是静态是初始化的定义方法：

```Java
//方法1
int[] intArr //定义好一个int数组类型的变量，也就是引用变量
intArr = new int[] {5,6,7,8} // 将这个引用变量指向一个int数组
// 方法2，两步一起来
int[] intArr = {5,6,7,8}
```

动态初始化方法

```Java
// 方法1
int[] prices
price = new int[5]
// 方法2
price = new int[5]
```

### 数组使用

数组的常见用法就是访问数组元素，包括对数组元素进行赋值和取出数组元素的值。可以单独访问一个元素，也可以通过for循环和foreach循环的方式对数组进行遍历。当然`foreach`是Java5之后出现的特性，举个例子

```Java
public class forEachTest
{
    public static void main(String[] args)
    {
        int[] prices;
        prices = new int[5];
        for (int i =0; i < prices.length; ++i)
        {
            prices[i] = i;
        }
        //price是形参，prices则是数组名
        for (int price : prices)
        {
            System.out.println(price);
        }
    }
}
```

**注**：foreach只能用于遍历数组，而不能对数组赋值。

### 数组的一些机制

为了更好理解的数组，你需要知道两个新的概念：堆(heap)和栈(stack)。

>**栈**（stack）又名堆栈，它是一种运算受限的线性表。其限制是仅允许在表的一端进行插入和删除运算。这一端被称为栈顶，相对地，把另一端称为栈底。向一个栈插入新元素又称作进栈、入栈或压栈，它是把新元素放到栈顶元素的上面，使之成为新的栈顶元素；从一个栈删除元素又称作出栈或退栈，它是把栈顶元素删除掉，使其相邻的元素成为新的栈顶元素。
>**堆**(heap):是计算机科学中一类特殊的数据结构的统称。堆通常是一个可以被看做一棵树的数组对象。

我觉得你肯定会一脸懵逼的，所以看下图直观感受下数组在内存的存放形式。

![](http://oex750gzt.bkt.clouddn.com/17-12-14/8597191.jpg)

数组分为数组引用变量和数组的实际对象，这两个存放在不同的位置。编程时通过数组引用变量对实际的数组对象进行修改。

基本数据类型的数组定义，初始化和数组操作的代码和内存中的变化如下

```Java
public class BasicArrayTest
{
    public static void main(String[] args)
    {
        int[] numArry;
        numArry = new int[8]
        for (int i = 0; i < numArry.length; ++i)
        {
            numArry[i] = i + 1;
            System.out.println(numArry[i]);
        }
    }
}

```

![](http://oex750gzt.bkt.clouddn.com/17-12-14/24922895.jpg)

### 数组的拓展和用途

Java8增加一个工具类Arrays用于处理数组，可以在<https://docs.oracle.com/javase/8/docs/api/index.html>学习用法。

数组的用途也很广，比如说你可以编一个命令行围棋工具,但是目前的用途就是做题目。在知乎上有一个提问是一行Python能实现什么丧性病狂的功能，[几个小例子告诉你, 一行Python代码能干哪些事](https://zhuanlan.zhihu.com/p/23321351). 我们要根据本次学习的Java基本语法实现

- 输出特定字符"Love"拼成的心形
- 输出Mandelbrot图像
- 打印九九乘法表
- 输出斐波那契数列

让我们用Java代码根据Mandelbrot的定义利用数组生成一副Mandelbrot图.先说一件很尴尬的事情，Java的基本数据类型里面没有复数，且我学Java时间合起来也不到3天，所以我就通过搜索引擎找到了一个复数类凑合用。代码如下，比起代码更重要的是对数学公式的理解。

```Java
public class Mandelbrot
{
    public static void main(String[] args)
    {
        // define the iteration
        int maxIterations = 100;
        // define the position range
        int xstart = -80;
        int xend   =  40;
        int ystart = -40;
        int yend   =  40;
        // define the scale
        double scalingFactor = 1.0 / 40;
        //define the two dimension array
        char[][] mandelbrot = new char[xend - xstart][yend-ystart];
        for (int x = xstart; x < xend; ++x)
        {
            for (int y = ystart; y < yend; ++y)
            {
                // whether the mod in the iteration will be large than 2
                Complex current = new Complex(x * scalingFactor, y * scalingFactor);
                Complex temp = current;
                for (int iter = 0; iter < maxIterations; ++iter)
                {
                    temp = temp.times(temp).plus(current);
                    if (temp.mod() > 2){
                        mandelbrot[x-xstart][y-ystart] = ' ';
                        break;
                    }else
                    {
                        mandelbrot[x-xstart][y-ystart] = '*';
                    }
                }
            }
        }
        // output the result
         for (int y = ystart; y < yend; ++y)
        {
            for (int x = xstart; x < xend; ++x)
            {
                System.out.print(mandelbrot[x-xstart][y-ystart]);
            }
            System.out.print('\n');
        }
    }
}
```

结果咋那么丑呢。。

![](http://oex750gzt.bkt.clouddn.com/17-12-14/79606310.jpg)

参考资料:

- [浅析Mandelbrot集合及其图形的绘制](https://www.cnblogs.com/anderslly/archive/2008/10/10/mandelbrot-set-by-fsharp.html)
- [Java 复数类](https://www.math.ksu.edu/~bennett/jomacg/c.html)