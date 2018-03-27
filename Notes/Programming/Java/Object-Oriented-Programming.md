---
title: Java学习笔记(三)面向对象编程-以构造复数类为例
tags: Java, OOP
notebook: Java学习笔记
---
# Java学习笔记(三)面向对象编程-以构造复数类为例

## 什么是面向对象

## 如何定义类

Java定义类的基本方法如下

```Java
[修饰符] class 类名
{
    (>=0)构造器定义
    (>=0)成员变量
    (>=0)方法
}
```

其中类名的修饰符为public, final, abstract中任何一个，不写也行。类内部构造需要**构造器定义**，**成员变量**，**方法**。当然不写也行，只不过没有意义而已。

- 成员变量定义: `[修饰符] 类型 成员变量名 [=默认值];`其中修饰符可以为：public, protected, private(前面三选一), static, final.
- 方法定义如下，修饰符可以是public, protected, private(前面三选一), static, final, abstract(后者二选一)

```Java
[修饰符] 方法返回值类型 方法名(形参列表)
{
    //可执行语句组成的方法体
}
```

- 构造器定义， 修饰符可以是public,protected, private。 但是构造器名必须和类名相同，并且没有返回类型。

```Java
[修饰符] 构造器名(形参列表)
{
    //可执行语句组成构造器执行体
}
```

以复数类构造为例说明下

```Java
public class Complex extends Object
{
    // 复数构造器， z = u + i * v
    public Complex(double u, double v)
    {
        x = u;
        y = v;
    }
    // 成员变量，两个double类型
    private double x,y;
    // 方法之一：返回实部
    public double real()
    {
        return x;
    }
    // 方法之一：返回虚部
    public double imag()
    {
        return y;
    }
}
```

这几行代码的做的事情如下：

- 首先定义复数构造器，构造器里面有两个形参。
- 随后定义了两个成员变量，再复数构造器中被形参赋值。由于类内部的成员无先后顺序，因此看起来好像没有定义成员变量就开始赋值，其实没有问题。
- 接着就是复数类的两个方法，返回实部和虚部

那么构造的复数类应该如何使用呢？其实和基本数据类型定义一样。

```Java
// 先定义后赋值
Comoplex cmpl;
cmpl = new Complex(1,2);
// 同时完成定义和赋值
Complex cmpl = new Complex(1,2);
```

将上述代码保存到Complex.java，然后新建一个ComplexTest.java， 测试之前构造的复数是否能用。

```Java
public class ComplexTest
{
    public static void main(String[] args)
    {
        Complex cmpl;
        cmpl = new Complex(1,2);
        System.out.format("real is %f, imag is %f", cmpl.real(), cmpl.imag());
    }
}
```

上面这部分内容如果要用C语言实现的话，我想到了结构体。

```C
#include <stdio.h>
struct Complex
{
    double real;
    double imag;
};
double Real(struct Complex c)
{
    return c.real
}
double Image(struct Complex c)
{
    return c.imag
}
int main(void)
{
    double Real(struct Complex c);
    double Image(struct Complex c);
    struct Complex cmpl = {1,2};
    printf("the image of comple 1 + 2i is %f", Real(cmpl));
}
```