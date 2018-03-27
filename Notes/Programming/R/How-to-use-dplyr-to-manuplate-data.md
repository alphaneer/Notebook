# 如何使用dplyr处理数据

以下操作基于R语言自带数据集"mtcars"和"iris".

## 列选取: select

让我们先从最简单的开始,`select`从数据框中提取其中几列。有两种方式，第一种是直接指定变量名

```r
# 选择hp列
select(mtcars, hp)
select(mtcars, HP=hp) #得到的列名为HP
# 选择hp到wt列
select(mtcars, hp:wt)
# 反向选择, 剔除hp的其他列
select(mtcars, -hp)
```

第二种是通过如下几个函数选择变量命中包含某个特征的列

- starts_with(): 以某个字符串开头
- ends_with(): 以某个字符串结尾
- contains(): 含有某个字符串
- matches(): 使用正则表达式去匹配
- num_range(): 选择数字结尾
- one_of(): 从提供的变量中选择
- everything(): 全选

比如说在iris数据集中选择以"Se"开头的列，或者以"th"结尾的列

```r
select(iris, starts_with("Se"))
select(iris, ends_with("th"))
```

以上基本能符合日常需求，当然你有一些奇怪的要求，比如说你想选择都小于4的列，可以用`select_if`

```r
select_if(mtcars, ~all(.<=4))
```

全选还能使用`select_all(mtcars)`，尽管这等价于`select(mtcars, everythins())`。

## 行选取: filter

dplyr提供`filter`,`filter_at`,`filter_all`,`filter_if`用来对行数据进行过滤，仅保留符合要求的行。

`filter`用来处理特定的几个变量，比如说:

```r
# 保留mtcars中cyl为6的行
filter(mtcars, cyl==6)
# 保留cyl为6但是vs不为0的行
filter(mtcars, cyl==6, vs != 0)
```

功能比较简单，适用于处理特定几列。但是如果你想要在包含10多个变量的数据框中找到包含NA的行，肯定不能逐个变量写判断语句，这就需要用到`filter_all`。

```r
filter_all(mtcars, any_vars(is.na(.)))
```

这里的`any_vars`表示任意一列，`.`则指代当前选定列，那么`is.na(.)`就会返回当前列的布尔值向量。

或者如果你的变量命命名很有特点，那么可以使用`filter_at`判断变量名是否符合要求来选择列：

```r
#判断变量名里有‘d’的列是否为NA.
filter_at(mtcars, vars(contains('d')), any_vars(is.na(.)))
# 判断变量名中含有'd'的列中，所有值是否都大于3
filter_at(mtcars, vars(contains('d')), all_vars(.>3))
```

这里的`vars()`函数用于选择变量，可以用`select`里用到筛选函数对变量名进行判断, 比如说`contains()`就返回变量命包含某个字符的变量。

更复杂一点，如果我们仅需要对那些整数列进行筛选，那么就需要用到`filter_if`, 因为它的第二个参数会传递给`rlang::as_function`用于构造闭包，而闭包则是一种函数。

```r
filter_if(mtcars, ~all(floor(.) == .), all_vars(.!=0))
```

在`~all(floor(.) ==.)`的参数选择下，最后是对cyl, hp, vs, am, gear, carb这几列进行处理，可以用如下的方式进行验证。

```r
filter_if(mtcars, ~ all(floor(.) == .), all_vars(. == 4))
```

由于hpb里不存在一个值为4，那么all_vars返回的是全FALSE的向量。