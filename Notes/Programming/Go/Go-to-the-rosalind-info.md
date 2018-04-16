# 用Go解决Rosalind练习题

## Introduction to the Bioinformatics Armory

统计核酸数目

```go
package main
import (
	"fmt"
	"io/ioutil"
	"log"
	"os"
)
func main() {
	args := os.Args
	fn := args[1:2]
	nucl, err := ioutil.ReadFile(fn[0])
	if err != nil {
		log.Fatal(err)
	}
	a_num := 0
	c_num := 0
	t_num := 0
	g_num := 0
	for i, n := 0, len(nucl); i < n; i++ {
		switch {
		case nucl[i] == 'A':
			a_num += 1
		case nucl[i] == 'G':
			g_num += 1
		case nucl[i] == 'T':
			t_num += 1
		case nucl[i] == 'C':
			c_num += 1
		default:
			log.Print(nucl[i])
		}
	}
	fmt.Printf("A\tC\tG\tT\n")
	fmt.Printf("%v\t%v\t%v\t%v\t\n", a_num, c_num, g_num, t_num)
}
```

## Introduction to Protein Databases

蛋白质数据库中心[UniProt](http://www.uniprot.org/)提供了蛋白详细的注释，如功能描述，功能与结构，翻译后修饰。它还支持蛋白相似性搜索，分类分析和文献引用等。

已知给定一个uniprot id,可以通过链接"<http://www.uniprot.org/uniprot/uniprot_id.txt>"或"<http://www.uniprot.org/uniprot/uniprot_id>" 获取关于该编号的详细描述。 通过编程的方式根据一个uniprot ID获取其参与的生物学进程(biological processes)

我使用Go的os.Args读取命令行参数中的编号，使用"net/http"获取响应, 使用"ioutil.ReadAll"获取响应中主体，返回字符数组。利用正则表达式进行解析，然后使用for循环提取出目标区段。

```go
package main

import (
    "net/http"
    "fmt"
    "os"
    "log"
    "io/ioutil"
    "regexp"
)

func main(){
    id := os.Args[1]
    link := "http://www.uniprot.org/uniprot/" + id + ".txt"
    resp, err := http.Get(link)
    if err != nil{
        log.Fatal(err)
    }
    content, err := ioutil.ReadAll(resp.Body)
    resp.Body.Close()
    if err != nil{
        log.Fatal(err)
    }
    re := regexp.MustCompile("P:(.*?);")
    BP := re.FindAllStringSubmatch(string(content[:]),-1)
    for i,n := 0, len(BP); i<n; i++{
        fmt.Println(BP[i][1])
    }
}
```

## GenBank Introduction

分子生物学家可获取的最大的整合型数据库就是[GenBank](https://www.ncbi.nlm.nih.gov/genbank/), 它包含了几乎所有公共的DNA序列和蛋白序列。 GenBank最早由NCBI在1982年建立，现在30多年过去了，里面存储的数据量超乎你的想象。

每个GenBank都有唯一的识别号，用于提取全序列，比如说`CAA79696,NP_778203, 263191547, BC043443, NM_002020`. 当然还可以用一些关键字搜索一类序列。

问题：GenBank包括如下几个子类，如Nucleotide, GSS(Genome Survey Sequence), EST(Expressed Sequene Tags). 为了精确从这些数据库中找到自己目标，需要用到一些搜索语法，比如说`(Drosophila[All Fields])`表示在所有区域中搜索Drosophila。 那么给定一个物种名，和两个日期，找到在这段时间内该物种上传到GenBank的氨基酸数。

解决方案：NCBI提供Entrez用于检索它存放的所有数据，并提供了相应的网页API<https://eutils.ncbi.nlm.nih.gov/entrez/eutils/>, 按照要求构建URL发起请求后就能返回目标响应用于解析。为了避免对服务器造成太大压力，NCBI对请求有一定的限制限制，每一秒不超过3个请求, 除非你在请求中带上了`API_KEY`。

> API key可以在<https://www.ncbi.nlm.nih.gov/account/>申请，允许每秒10个请求。考虑到国内这个网速，我觉得应该是用不到的
>![数据库](http://oex750gzt.bkt.clouddn.com/18-3-29/2662035.jpg)
> 可供查询的数据库

Entrez搜索语法：

- Boolean操作: AND OR NOT
- 限定领域: `[]`, 如 `horse[Organism]`
- 日期或其他范围: `:`， 如日期 `2015/3/1:2016/4/30[Publication Date]`, 如序列长度`110:500[Sequence Length]`

构建URL: 以样本数据`Anthoxanthum 2003/7/25 2005/12/27`为例，搜索语法为`Anthoxanthum[Organsim] AND 2003/7/25:2005/12/27[Publication Date]，转换成URL link就是

```bash
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=Anthoxanthum[Organsim]+AND+2003/7/25:2005/12/27[Publication Date]
```

> 如果搜索语法中有`"`和"#", 需要转换成"%22","%23", `[]`外的空格要用"+"代替,`[]`和`()`内空格要用"%20"替换。

解决这个问题不能想的太复杂，我们需要假设输入时是`"Anthoxanthum[Organsim]+AND+2003/7/25:2005/12/27[Publication Date]"`，而不是`Anthoxanthum 2003/7/25 2005/12/27`, 这样子我们就只需要构建URL link, 然后发送请求解析响应的xml就行

```go
//GBK.go
package main

import (
        "fmt"
        "io/ioutil"
        "log"
        "net/http"
        "os"
        "encoding/xml"
)

// the base url link of entrez API
const BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

func buildLink(db, query string) string {
        link := BASE + "?db=" + db + "&term=" + query
        return link
}

func eSearch(link string) {
        resp, err := http.Get(link)
        if err != nil {
                log.Fatal(err)
        }
        body, err := ioutil.ReadAll(resp.Body)
        resp.Body.Close()
        if err != nil {
                log.Fatal(err)
        }
}
type eSearchResult struct {
        Count  string   `xml:"Count"`
        RetMax string   `xml:"RetMax"`
        IdList []string `xml:"IdList>Id"`
}
func main() {
        db := os.Args[1]
        query := os.Args[2]
        link := buildLink(db, query)
        fmt.Println(link)
        content := eSearch(link)
        result := eSearchResult{}
        err := xml.Unmarshal(content, &result)
        if err != nil {
                fmt.Printf("err: %v", err)
                return
        }
        fmt.Printf("Count:%s\n", result.Count)
        fmt.Printf("Idlist:%v\n", result.IdList)
}
}
```

读取命令行的字符串，第一个为数据库，第二个为请求。 将参数传入后构建成link，使用"net/http"发起请求,将得到的响应用`ioutil.ReadAll`读取保存为**字符数组**. 然后定义结构体用来保存XML的解析结果。xml.Unmarshal的使用涉及到**interface**的知识。

## Data Formats

表示序列最常见的格式就是FASTA(.fas,. fasta), 当然还有NEXUS(.nex, .nexus, nxs)和PHYLIP(.phy)等，每一个格式都觉得自己挺好，而你又得用他们软件，所以很多时候就得面对一个问题，如何把A格式转成B格式。

GenBank是世界上最大的生物序列数据存储中心, 他们有一套自己的格式叫做GenBank, 用于记录这段序列的方方面, 适合人类阅读, 就是不太适合直接用于数据分析, 往往需要进一步转成FASTA才行。

问题：给定不多与10个的GenBank登记号，找到其中最短的序列，以FASTA的格式展示。

解题: 根据提供PubMed的entry使用esearch获取UID, 如下代码使用了之前写的程序'GBK'

```bash
cat rosalind_frmt.txt | tr ' ' '\n' | xargs -i ./GBK nuccore {} | sed -n '/Idlist/p' | cut -d '[' -f 2 | cut -d ']' -f 1 > UID.list | tr '\n' ' ' > uid.list
```

使用entrez的efetch功能下载序列，efetch可以指定format为fasta，从而下载fasta格式文本，下载之后涉及到IO操作，保存文件。

```go
package main
import (
	"fmt"
	"io/ioutil"
	"log"
	"net/http"
	"os"
	"strings"
)
// the base url link of entrez API
const BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
func buildLink(db, format string, query []string) string {
	s := strings.Join(query,",")
	link := BASE + "?db=" + db + "&id=" + s + "&rettype=" +format
	return link
}
func eFetch(link string) []byte {
	resp, err := http.Get(link)
	if err != nil {
		log.Fatal(err)
	}
	body, err := ioutil.ReadAll(resp.Body)
	resp.Body.Close()
	if err != nil {
		log.Fatal(err)
	}
	return body
}
func main() {
	db := os.Args[1]
	format := os.Args[2]
	query := os.Args[3:]
	link := buildLink(db,format, query)
	fmt.Println(link)
	content := eFetch(link)
	ioutil.WriteFile("sequence."+format, content,0660)
}
```

通过`./FRMT nuccore fasta id1 id2`下载序列，下载的序列用来比较长度。