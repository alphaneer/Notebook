# Pysam学习笔记

## 操作VCF/BCF文件

读取和写出

```Python
from pysam import VariantFile
bcf_in  = VariantFile("test_in.vcf", "r")
bcf_out = VariantFile("test_out.vcf", "w", header=bcf_in.header)
for rec in bcf_in.fecth():
    bcf_out.write(rec)
```

VariantFile函数得到的是 **pysam.libcbcf.VariantFile** 对象, 这是一个可遍历对象, 通过`dir()`可以发现它有`__iter__`和`__next__`方法。因此如果仅仅是遍历全部记录，那么`__iter__`等价于`fecth`.

```Python
type(bcf_in) # 对象类型
dir(bcf_out) # 方法
```

VCF格式分为Header和Record两个部分. record记录每个变异位点的具体信息，为了从中提取所需数据，需要理解Pysam的解析策略。

```Python
rec1 = bcf_in.__next__()
dir(rec1)
```

vcf的record每一行都是9列+N列样本(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample1, sample2,..), 解析之后就是如下方法

- .chrom: 返回字符串
- .pos: 返回数值。 这个是以0为基, 可以用.start和.stop
- .id: 如果无记录, 就是NoneType
- .ref: 返回字符串
- .alt: 返回元祖(tuple), 因为一个位点上可以有多个变异类型
- .qual: 返回数值
- .filter: 返回pysam.libcbcf.VariantRecordFilter对象, 类似于字典
- .info: 返回pysam.libcbcf.VariantRecordInfo对象，类似于字典, 存放所有样本的统计信息
- .format: 返回pysam.libcbcf.VariantRecordFormat，类似于字典, 存放后续每个样本数据存放顺序和数据类型
- .samples: 返回pysam.libcbcf.VariantRecordSamples, 类似于字典, 存放每一个样本的具体信息

.filter, .info, .format, .samples虽然都能返回类字典(或者说哈希表)数据结果，但是在方法上存在差别。

VariantRecordFilter对象可以通过`.filter.add`增加过滤类型, 当然需要事先在header中添加元信息，如下:

```Python
bcf_in.header.filters.add(id="ugly",number=None, type=None,description="i don't likt it") #增加员信息
rec = bcf_in.__next__()
rec.filter.add("ugly") # 增加过滤条件
rec.filter.keys() # 查看
```

VariantRecordInfo对象可以删除一个键值对(pop),可以更新已有的键值对。

```Python
rec.info.pop('TYPE') # 删除TYPE
rec.info['ODDS'] # 变更前
rec.info.update({'ODDS':12}) #变更
rec.info['ODDS'] # 变更后
```

VariantRecordFormat和VariantRecordSamples关系比较紧密，但前者只能查看不提供方法进行修改, 而VariantRecordSamples和VariantRecordInfo一致。由于可以有多个样本，提取数据的时候就需要多层迭代，例如提取所有样本的GT

```Python
for key,value in rec.samples.iteritems():
    print(key, value['GT'])
```

例如只有两个样本，我想比较这两个样本的GT是否相同

```Python
GT = [value['GT'] for value in rec.samples.values()]
GT[0].__eq__(GT[-1])
```

综上，就可以在Python中写出一个过滤器剔除缺失基因组记录，保留其中样本基因组纯合但不同的记录

```Python
import sys
from pysam import VariantFile as vcf

if len(sys.argv) < 3:
    sys.exit(1)
else:
    in_name  = sys.argv[1]
    out_name = sys.argv[2]

bcf_in  = vcf(in_name)
# add metadata
command = "##pysamCommand=GT[0].__ne__((None,)) and GT[-1].__ne__((None,)) and GT[0].__ne__(GT[-1]) and GT[0].__ne__((0,1)) and GT[-1].__ne__((0,1))"
bcf_in.header.add_line(command)
bcf_out = vcf(out_name, "w", header=bcf_in.header)

for rec in bcf_in.__iter__():
    GT = [value['GT'] for value in rec.samples.values()]
    if GT[0].__ne__((None,)) and GT[-1].__ne__((None,)) and \
           GT[0].__ne__((0,1)) and GT[-1].__ne__((0,1)) and \
           GT[0].__ne__(GT[-1]):
        bcf_out.write(rec)

```