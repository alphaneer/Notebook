# 学会使用Go标准库

## 以CSV库为例学习Go的标准库用法

标准库提供了一系列实用函数的API，让我们避免重复造轮子. 为了学会使用一个已有的轮子，我们需要学看文档。Go所有标准库文档都放在<https://golang.org/pkg/>, 假如我们想读取一个CSV文件, 我们就需要学习csv包的用法。

首先是**加载**， 即`import "encondig/csv". 这一步声明了csv源代码所在路径，即`/path/to/go/src/encondig/csv"

其次简单了解这个包需要提供什么样的输入，即参数，会返回什么样的输出，即返回的数据类型。CSV读写符合"RFC 4180"标准的要求，类似于"field1,field2,field3", 所以你要保证field中不应该存在额外的逗号,为了已经有了逗号，那么就需要在两端加上引号。对于被引号包围的字符串，会默认移除两边的引号，如果你想保留引号，就得用两个引号

```bash
## 输入
"the ""word"" is true","a ""quoted-field"""
"Multi-line
field","comma is ,"
## 输出
{`the "word" is true`, `a "quoted-field"`}
{`Multi-line
field`, `comma is ,`}
```

在换行符移除前进行归位(carriage return), 即回到一行字的开头，但是不代表换行。

> 关于换行符的故事，参考维基百科<https://en.wikipedia.org/wiki/Enter_key>

再接着会说明这个包定义了哪些变量，数据类型，以及返回这类数据的函数，并且每一个函数都有案例帮助理解。

```go
package main

import (
	"encoding/csv"
	"fmt"
	"io"
	"log"
	"strings"
)
// 导入多个标准包，fmt用于标准化输出,io为输入输出相关，log用于记录日志，strings是强化的字符类型
func main() {
	in := `first_name,last_name,username
"Rob","Pike",rob
Ken,Thompson,ken
"Robert","Griesemer","gri"
`
// 用``可以多行输入
	r := csv.NewReader(strings.NewReader(in))
//调用strings.NewReader函数读取字符串，返回Reader结构体的指针，而Reader实现io.Reader
//使用csv.NewReader函数读取io.Reader, 返回csv包定义的Reader结构体的指针
	for {
        record, err := r.Read()
        //Reader实例的Read方法
        //一般情况下返回字符数组和ErrFieldCount(用于判断每一行的列数是否都相同)
        //否则返回的err为non-nil record或non-nil error
        //没有数据返回io.EOF
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatal(err)
		}

		fmt.Println(record)
	}
}
```

案例为了方便说明，于是将要读取的csv格式的字符串放在了代码中，实际情况下可以将这部分代码保存为'test.csv'，然后想办法给 strings.NewReader 提供字符串数据类型。我的策略是用`ioutil.ReadFile`读取所有字符得到字符数组`[]byte`，然后将字符数组显式转换字符串数据类型。

```go
package main

import (
	"encoding/csv"
	"fmt"
	"io/ioutil"
	"log"
	"strings"
	"io"
)

func main() {
	dat, err := ioutil.ReadFile("test.csv")
	if err != nil {
		log.Fatal(err)
	}
	r := csv.NewReader(strings.NewReader(string(dat[:])))
	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatal(err)
		}
		fmt.Println(record)
	}
}
```

## 常用标准库

命令行参数: flag, os.Args