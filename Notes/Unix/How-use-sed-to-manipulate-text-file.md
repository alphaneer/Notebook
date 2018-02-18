# 使用sed处理文本

sed语法

```bash
sed [options] {sed-commands} {input-file}

sed [options] -f {sed-commands-in-a-file} {input-file}

sed [options] -e {sed-command-1} -e {sed-command-2} {input-file}

sed
 [options] '{
sed-command-1
sed-command-2
}' input-file

```

SED流： 读 -> 执行 -> 打印 -> 重复

源文件

```bash
$ cat source.txt
101,Ian Bicking,Mozilla
102,Hakim El Hattab,Whim
103,Paul Irish,Google
104,Addy Osmani,Google
105,Chris Wanstrath,Github
106,Mattt Thompson,Heroku
107,Ask Solem Hoel,VMware
```

范围

```bash
sed -n '2~3 p' source.txt
```

模式匹配

```bash
sed -n '/Paul/ p' source.txt
sed -n '/Paul/,5 p' source.txt
sed -n '/Paul/,/Addy/ p' source.txt
sed -n '/Paul/,+2 p' source.txt
```

删除行

```bash
sed '2 d' source.txt
sed '/^$/ d' source.txt
sed '/^#/ d' source.txt
```

重定向

```bash
sed 'w output.txt' source.txt
sed -n '1,4 w output.txt'  source.txt
sed -n '/Ask/,$ w output.txt'  source.txt
```

替换

```bash
sed '[address-range|pattern-range] s/original-
string/replacement-string/[substitute-flags]' inputfile
```

行后增加

```bash
sed '[address] a the-line-to-append' input-file
```

行前增加

```bash
sed '[address] i the-line-to-insert' input-file
```

高级话题：

- 处理多行模式空间（N、D、P）。
- 采用保持空间来保存模式空间的内容并使它可用于后续的命令（H、h、G、g、x）。
- 编写使用分支和条件指令的脚本来更改控制流（：、b、t）。

前面基本用法中也有提到模式空间，即为处理文件中一行内容的一个临时缓冲区。处理完一行之后就会把模式空间中的内容打印到标准输出，然后自动清空缓存。

而这里说的保持空间是sed中的另外一个缓冲区，此缓冲区正如其名，不会自动清空，但也不会主动把此缓冲区中的内容打印到标准输出中。而是需要以下sed命令进行处理：

      d     Delete pattern space.  Start next cycle.    删除pattern space的内容，开始下一个循环.
      h、 H    Copy/append pattern space to hold space.   复制/追加pattern space的内容到hold space.
      g、 G    Copy/append hold space to pattern space.   复制/追加hold space的内容到pattern space.
      x      Exchange the contents of the hold and pattern spaces.    交换hold space和pattern space的内容.

正元字符集

```bash
^ 匹配行开始，如：/^sed/匹配所有以sed开头的行。
$ 匹配行结束，如：/sed$/匹配所有以sed结尾的行。
. 匹配一个非换行符的任意字符，如：/s.d/匹配s后接一个任意字符，最后是d。
* 匹配0个或多个字符，如：/*sed/匹配所有模板是一个或多个空格后紧跟sed的行。
[] 匹配一个指定范围内的字符，如/[ss]ed/匹配sed和Sed。
[^] 匹配一个不在指定范围内的字符，如：/[^A-RT-Z]ed/匹配不包含A-R和T-Z的一个字母开头，紧跟ed的行。
\(..\) 匹配子串，保存匹配的字符，如s/\(love\)able/\1rs，loveable被替换成lovers。
& 保存搜索字符用来替换其他字符，如s/love/**&**/，love这成**love**。
\< 匹配单词的开始，如:/\<love/匹配包含以love开头的单词的行。
\> 匹配单词的结束，如/love\>/匹配包含以love结尾的单词的行。
x\{m\} 重复字符x，m次，如：/0\{5\}/匹配包含5个0的行。
x\{m,\} 重复字符x，至少m次，如：/0\{5,\}/匹配至少有5个0的行。
x\{m,n\} 重复字符x，至少m次，不多于n次，如：/0\{5,10\}/匹配5~10个0的行。
```