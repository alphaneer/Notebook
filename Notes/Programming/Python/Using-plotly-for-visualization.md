# 使用plotly进行数据可视化

## plotly起步

```bash
pip install plotly
# 升级
pip install plotly --upgrade
```

plotly作图结果有两种展示方式，一种是在线托管作图（仅当私人托管时才需要付费）；另一种是生成一个离线文件用于查看。后者所需要的函数是`plotly.offline.plot()`和`plotly.offline.iplot()`，分别生成一个正常的网页文件和或在Jupyter notebook里展示。

```Python
# html
import plotly
from plotly.offline import plot
import plotly.graph_objs as go

plot({
    "data": [go.Scatter(x=[1,2,3,4],y=[4,3,2,1])],
    "layout": go.Layout(title="Hello World!")
})
# Jupyte Notebook下
import plotly
from plotly.offline import iplot
import plotly.graph_objs as go
plotly.offline.init_notebook_mode(connected=True)
iplot({
    "data": [go.Scatter(x=[1, 2, 3, 4], y=[4, 3, 2, 1])],
    "layout": go.Layout(title="hello world")
})
```

由于这两种方法接受相同的plotly.graph_objs数据类型，所以就可以先在Jupyter Notebook先进行数据探索整理代码，后续就能自动化处理。

## 用plotly做出能自定义控制的图