# Go的基本语言特性

## 数据类型

数据类型分为两类: 基本类型和引用类型。基本数据类型包括布尔类型，数值类型、浮点类型和字符类型等。引用类型包括slice, map和channel，有着复杂的内存结构，需要申请内存以及初始化相关属性。

```go
// 基本数据类型
var x,y,z int
var f = 1.6
// 引用类型
a := []int{0,0,0} // 提供初始化表达式
b := make([]int,3) //slice
```

## 函数

Go的函数不支持嵌套(nested)、重载(overload)和默认参数(default parameter), 不过在Go的函数是第一类对象，可以作为参数进行传递。

## 数据结构

### 数组

```go
a := [3]int{1,2} // 未初始化的元素为0
c := [5]int{2:100, 4:200} // 初始化部分值
```

内置函数`len`和`cap`返回数组长度

### 结构体

C语言只有数据类型就是数组，其他数据类型就靠结构体进行定义。Go的结构体和C的差不多，支持自身类型的指针成员，也就是说它能搞出链表了。

```go
package main

type entry struct {
	value int
	next  *entry
}

func main() {
	n1 := entry{value: 100}
	n2 := entry{value: 200, next: &n1}
	n3 := entry{value: 300, next: &n2}
	println(n3.next.next.value)

}
// 猜猜看，最后输出结果是啥
```

### 切片

切片通过内部指针和相关属性引用数组片段，从而实现变长。

```go
data := [...]int{0,1,2,3,4,5,6}
slice := data[1:4:5] //[low:high:max]
```

对切片对象使用`len`和`cap`和定义时的low,high,max有关。len=high-low, cap=max-low. len表示从数组获取的元素数量，切片的读写操作只能在这个范围进行，而cap则是切片可以变化的范围，但是不能超过原数组大小。

可以不通过数组，直接创建切片对象

```go
s1 := make([]int, 6, 8)
fmt.Println(s1, len(s1), cap(s1)
```

### Map

Map是一种引用类型，类似于Python的字典, R语言的列表，其实都是哈希表(hash table)。Map是一种键值对关系，其中键必须是支持相等运算符号的数据类型，而值可以是任意类型。几个定义Map的例子

```go
m := map[string]int{
    "a":1
}
// 这里的key是字符, value是整数值
m := map[int]struct{
    name string
    age int
}{
    1: {"user1",10},
    2: {"user2",20}
}
//这里的key是整数值, value是结构体，看起来比较复杂，可以分为两个部分
type user struct{
    name string
    age int
}
m := map[int]user{
    1: {"user1",10},
    2: {"user2",20}
}
```

## 结构体与面向对象

了解了结构体很容易让我联想到面向对象编程，Go语言不算是完全的OOP语言，它仅仅支持OOP三大特征里的封装和通过匿名字段实现的类继承，而没有多态。

```go
type User struct{
    id int
    name string
}

func (self *User) TestPointer(){
    fmt.Printf("TestPointer: %p, %v\b", self, self)
}
func (self User) TestValue(){
    fmt.Printf("TestValue: %p, %v\n", &self, self)
}
func main(){
    u := User{1, "Tom"}
    fmt.Printf("User: %p, %v", &u,u)
    mv := User.TestValue
    mv(u)

    mp := (*User).TestPointer
    mp(&u)

    mp2 := (*User).TestValue
    mp2(&u)
}

```
