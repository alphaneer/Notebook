# Python的正则表达式

Python通过导入标准库re实现正则表达式(regular expression)，Python的正则表达式引擎和Perl一样，并且兼容Perl流派的元字符。

## 元字符

Python支持的元字符很多，一种是比较常见，我之前也就只会用这些

- `.`表示任意一个字符，默认不匹配换行符，制表符
- `|`表示或，`ca|bd`会匹配ca或bd,而不是cab, cbd, 如果想要匹配后者，则需要用到`()`进行分组
- `^`,`$`表示位置符号，行首和行尾 如`^ab$`匹配ab, 不匹配eab, abe,aeb
- 量词，表示重复数，`*`任意多次, `+`一次以上, `?`0次或一次, `{m,n}`m~n次, `{m}`重复m次,`{m,}`重复大于m次
- 在上述量词后接`?`, 就从贪婪模式变为非贪婪模式。举个例子，对于`abbbbbb`这个字符串，`ab*`和`ab*?`的结果不同，前者匹配`abbbbbb`，后者匹配`a`，也就是贪婪模式尽可能多匹配。
- `[...]`表示多选项，比如`a[bc]`就可以匹配ab,ac, 如果是`[a-z]`那么表示从a到z范围. 所有元字符在`[]`中都会被认为是普通字符。所有元字符在`[]`
- `(...)`表示捕获型分组，被`(...)`匹配到部分，可以用`\1`,`\2`进行引用
- "\" 表示转义，由于该符号也是字符串的元字符，那么在构建模式的时候要万分小心，因为Python会先对字符串进行加工，然后才会传入到正则引擎中。也就是说，也就是如果你想匹配"\" , 你的模式写法得是`\\\\`，因为如果只写`\\`,会被Python先翻译成`\`,所以必须写成`\\\\`。因此建议用使用原始字符串(raw string),即`r"\\"`

下面的一些比较高级，在我写作时能记得的元字符，基本上都是`(?...)`一类的增强型标记，具体含义和`?`后紧接的第一个字符有关

- `(?:...)`: 非捕获型分组，也就是仅仅分组，正则引擎不会记住他用于后续引用
- `(?=...)`: 向后检查，要求当前位置后符合`...`表示的模式, `(?!...)`也是向后检查，只不过要求当前位置紧接的内容不能被`...`匹配
- `(?<=...)`和`(?<!...)`是向前检查。

在《精通正则表达式》中，作者举了一个例子，将"12345679"变为更容易阅读的"12,345,679"形式。 也就是找到一个位置前面是数字，后面是3的倍数个数字的位置插入逗号

```Python
re.sub(r"(?<=\d)(?=(\d\d\d)+$)",",","1234567")
```

下面是我需要翻阅资料才能记得

- `(?P<name>...)`: 在之前捕获型括号的基础上，将捕获到的内容赋值给`name`, 其中该内容可以用`(?P=name)`进行引用
- `(?#...)`: 这个仅仅是注释，不做任意匹配
- `(?aiLmsux)`比较复杂，记不太起来
- `(?(id/name)yes-pattern|no-pattern)`更加复杂，需要举一个例子。`(<)?(\w+@\w+(?:\.\w+)+)(?(1)>|$)`来解释，当然这个例子理解起来也不容易。解释起来就是，第一个括号先尝试**捕获**匹配`<`, 编号为1，然后是第二个括号匹配“字符串@字符串”，比如说user@host,然后第三个括号表示不捕获分组, 识别".com"这类，然后第四个括号就是看第一个括号有没有捕获到东西，如果有就去匹配`>`，没有则是匹配行尾。也就是你的邮箱地址要么为"user@host.com",要么为`<user@host.com>`,其他都是不符合要求。

## 常用函数

一般用法都是用`re.compile`构建一个正则表达式对象，这个正则表达式对象可以用在`re.match`,`re.search`,`re.find`,`re.findall`等函数里，同时该对象也有`.match`,`.search`方法。举个例子，比如说你知道了一个形如GSExxx的GEO编号，你需要提取这个编号下的所有GSMxxx编号，然后根据这个GSMxxx编号去提取SRA编号，以随便找的GSE100566为例。

首先利用Python的requests库抓取网页信息

```Python
# Python
import re
import requests
base_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
acc      = "GSE101571"

resp = requests.get(base_url + acc)
contents = resp.text()
```

然后构建一个正则表达式，去捕获所有的GSMxxx类型的编号

```bash
pattern = "GSM\d+"
GSM_acc = re.findall(pattern, contents)
```

或许你不满足于此，你还希望捕获到每个GSM编号后的描述，也就是"GSM2686880	SET-2_STAT1-D1",这两个部分你都需要。通过检查网页元素，你发现了一个规律，也就是这两个元素是在一个`tr`内

![HTML结构](http://oex750gzt.bkt.clouddn.com/18-5-8/30817865.jpg)

你信心慢慢的构建了一个匹配模式，结果啥都没有匹配到

```bash
pattern = re.compile("<tr><td.*?><a.*?>(GSM\d+)</a></td><td.*?>(.*?)</td>")
re.search(pattern, contents)
```

你发现这似乎由于这个HTML里有很多神奇的空白和"\n",原本方便人类阅读的记号却阻碍了数据处理，你必须做点什么，你想到了可以用`re.sub`进行替换，所以你做了如下的事情

```Python
contents = re.sub(r"\n\s*","",contents)
```

最后你终于用原来的匹配模式得到了以元组数据结构的结果

```Python
result = re.findall(pattern, contents)
```

下一步根据GSMxxx编号去提取SRX编码。这一步的核心就是从元祖中提取元素，然后构建一个url去爬取新的网页，然后提取SRX编号即可以。先测试第一个，

```Python
r1 = results[0][0]
r1_resp = requests.get(base_url + r1)
m = re.search("SRX\d+", r1_resp.text)
m.group(0)
```

然后开始遍历,存储到字典中。考虑到网络延迟所耽误的时间远远大于内存分配的时间，也就没有必须要预先分配内存空间。

```Python
sra_dict = {}
for acc in results:
    key = acc[0]
    resp = requests.get(base_url + key)
    value = re.search("SRX\d+",resp.text).group(0)
    sra_dict[key] = value
```

使用`re.match`, `re.search`, `re.fullmatch`和`re.finditer`都会返回匹配对象(match object).

```Python
match = re.search(pattern, string)
```

我们可以用`match.group()`从中提取`()`捕获到分组, 默认是0，也就是所有匹配分组，也可以以切片的形式，即`match[0]`,因为它有`.__getitem__`方法。如果用`mathc.groups()`,则范围是元组。

此外，还可以通过`match.start()`和`match.end()`获取匹配起始和结束为止

## 例子

正则可以方便的用在解析器中，比如说你官方的标准库就给了一个例子，关于语法解析。

```Python
import collection
import re

Token = collection.namedtuple('Token', ['typ','value','line','column'])

def tokenize(code):
    keywords = {'IF', 'THEN','ENDIF','FOR','NEXT','GOSUB','RETURN'}
    token_specification = [
        ('NUMBER', r'\d+(\.\d*)?'), # 整数或小数
        ('ASSIGN', r':='),          # 赋值负号
        ('END', r';'),              # 表示式结束
        ('ID', r'[A-Za-z]+'),       # 变量名
        ('OP', r'[+\-*/]'),         # 数学运算
        ('NEWLINE',r'\n'),          # 换行符
        ('SKIP', r'[ \t]+'),        # 跳过空白符号
        ('MISMATCH', r'.'),         # 其他不符合要求的符号
    ]
    tok_regex = '|'.join('(?P<%s>%s)' % pair for pair in token_specification)
    line_num = 1
    line_start = 0
    for mo in re.finditer(tok_regex,code):
        kind = mo.lastgroup
        value = mo.group(kind)
        if kind == 'NEWLINE':
            line_start = mo.end()
            line_num +=1
        elif kind == 'SKIP':
            pass
        elif kind == 'MISMATCH':
            raise RuntimeError(f'{value!r} unexpected on line {line_num}')
        else:
            if kind == "ID" and value in keywords:
                kind = value
            column = mo.start() - line_start
            yield Token(kind, value, line_num, column)

statements = '''
    IF quantity THEN
        total := total + price * quantity;
        tax := price * 0.05;
    ENDIF;
'''

for token in tokenize(statements):
    print(token)

```