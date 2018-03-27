# 使用dash开发交互式网页

后续的操作前,需要安装如下Python包

```bash
pip install dash==0.20.0  # The core dash backend
pip install dash-renderer==0.11.2  # The dash front-end
pip install dash-html-components==0.8.0  # HTML components
pip install dash-core-components==0.18.1  # Supercharged components
pip install plotly --upgrade  # Plotly graphing library used in examples
```

## Dash应用布局

### 使用Dash生成HTML

Dash应用包括两个部分，应用布局(layout)和数据交互(interactivity)。其中布局部分用来展示数据以及引导使用者使用。Dash提供了`dash_core_components`和`dash_html_components`, 以类的方式对HTML和JS进行封装，便于调用。下面先构建一个最简单的布局

```Python
import dash
import dash_core_components as dcc
import dash_html_components as html

app = dash.Dash()

app.layout = html.Div(children=[
    html.H1(children = 'Hello Dash'),
    html.Div(children = '''
        Dash: A web application frameworkd for Python.
        '''),
    dcc.Graph(
        id = 'example-graph',
        figure = {
            'dash':[
                {'x': [1,2,3], 'y':[4,1,2], 'type':'bar', 'name':'SF'},
                {'x': [1,2,3], 'y':[2,4,5], 'type':'bar', 'name':'Montrel'},
                ],
            'layout':{
                'title':'Dash data Visualization'
                }
            }
        )
])

if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0')
```

首先用`app=dash.Dash()`创建了Dash应用的实例，这个实例可以通过`app.run_server()`运行。

其次这个应用的布局(layout)由html组件(html.Div等)和图形组件(dcc.Graph等)构成。其中基础的html标签来自于`dash_html_components`，而更加React.js库里的高级组件则是由`dash_core_components`提供。

最后的展示形式需要后期慢慢的调整， 比如说调整一下字体对齐, 字体颜色和背景颜色等

```Python
import dash_core_components as dcc
import dash_html_components as html

app = dash.Dash()

colors = {
        'background':'#111111',
        'text':'#7FDBFF'
}

app.layout = html.Div(style={'backgroundColor':colors['background']},
    children=[
    html.H1(
        children = 'Hello Dash',
        style = {
            'textAlign':'center',
            'color': colors['text']
            }
        ),
    html.Div(children = '''
        Dash: A web application frameworkd for Python.
        ''', style = {
            'textAlign':'center',
            'color': colors['text']
            }
        ),
    dcc.Graph(
        id = 'example-graph',
        figure = {
            'data':[
                {'x': [1,2,3], 'y':[4,1,2], 'type':'bar', 'name':'SF'},
                {'x': [1,2,3], 'y':[2,4,5], 'type':'bar', 'name':'Montreal'},
                ],
            'layout':{
                'plot_bgcolor': colors['background'],
                'paper_bgcolor': colors['background'],
                'font':{
                    'color': colors['text']
                    },
                'title':'Dash data Visualization'
                }
            }
        )
])

if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0')
```

这里的html组件都设置了`style`，用来调整样式，

### 可视化

`dash_core_components`库中有一个`Graph`组件，它利用开源的JavaScript图形库--`plotly.js`进行交互式数据渲染。Graph里的`figure`参数等价于`plotly.py`里的`figure`参数，即任何`plotly.js`支持的图形都可以用`dash_core_components`调用。查看<https://plot.ly/python/>了解更多`plotly.py`的图形。

比如说这里可以基于Pandas的数据库创建散点图

```Python
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import pandas as pd

app = dash.Dash()

df = pd.read_csv(
    'https://gist.githubusercontent.com/chriddyp/' +
    '5d1ea79569ed194d432e56108a04d188/raw/' +
    'a9f9e8076b837d541398e999dcbac2b2826a81f8/'+
    'gdp-life-exp-2007.csv')

plot = [dcc.Graph(
        id = 'life-exp-vs-GDP',
        figure = {
            'data':[
                go.Scatter(
                    x=df[df['continent'] == i]['gdp per capita'],
                    y=df[df['continent'] == i]['life expectancy'],
                    text=df[df['continent'] == i]['country'],
                    mode='markers',
                    opacity=0.7,
                    marker={
                        'size':15,
                        'line':{'width':0.5, 'color':'white'}
                    },
                    name = i
                ) for i in df.continent.unique()
            ],
            'layout': go.Layout(
                xaxis={'type':'log','title':'GDP per Capita'},
                yaxis={'title':'Life Expectancy'},
                margin={'l':40,'b':40,'t':10,'r':10},
                legend={'x':0, 'y':1},
                hovermode='closest'
            )
        }
    )]

app.layout = html.Div(
    html.Div(children=[
        html.Div(className='col-md-4'),
        html.Div(plot,className='col-md-4')],
        className='row'
    )
)

# Append an externally hosted CSS stylesheet
my_css_url = "https://cdn.bootcss.com/bootstrap/3.3.7/css/bootstrap.min.css"
app.css.append_css({
    "external_url": my_css_url
})

# Append an externally hosted JS bundle
my_js_url = 'https://cdn.bootcss.com/bootstrap/3.3.7/js/bootstrap.min.js'
app.scripts.append_script({
    "external_url": my_js_url
})

if __name__ == '__main__':
    app.run_server(debug=True)
```

这部分代码将图形部分的代码从html组件中抽离出来，写完之后，再添加到html总体组件中。此外还增加了`bootstrap`的css样式。

### Markdown语法

Dash的`dash_html_components`支持原生的HTML语句，而`dash_core_components`得`Markdown`提供了Markdown得渲染。

```Python
import dash
import dash_core_components as dcc
import dash_html_components as html

app = dash.Dash()

markdown_text = '''
### Dash and Markdown

Dash apps can be written in Markdown.
Dash uses the [CommonMark](http://commonmark.org/)
specification of Markdown.
Check out their [60 Second Markdown Tutorial](http://commonmark.org/help/)
if this is your first introduction to Markdown!
'''

app.layout = html.Div([
    dcc.Markdown(children=markdown_text)
])

if __name__ == '__main__':
    app.run_server(debug=True)
```

`dash_core_components`里不仅仅提供了Markdown, graphs这些图形组件，还支持下拉栏等其他使用工具，可在<https://plot.ly/dash/dash-core-components>进一步了解

### 第一部分小结

这部分主要是学习了Dash应用得`layout`. `layout`是不同组件的层级关系树，最后结果是html页面。html页面的HTML基本语法由`dash_html_components`提供，而高级的绘图和下拉栏等则是由`dash_core_components`提供.

参考资料：

- <https://plot.ly/dash/getting-started>
- <https://plot.ly/dash/dash-core-components>
- <https://plot.ly/dash/dash-html-components>

## 响应式编程

第一部分完成了整体布局，但是基本都是静态图形，无法体现dash交互性数据探索特性。这一部分则是让图形能够动起来，对我们的操作有所回应。

```Python
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html

app = dash.Dash()

app.layout = html.Div([
    dcc.Input(id='my-id', value='initial vale', type='text'),
    html.Div(id='my-div')
])

@app.callback(
    Output(component_id='my-div', component_property='children'),
    [Input(component_id='my-id', component_property='value')]
)
def update_output_div(input_value):
    return 'you\'ve entered "{}"'.format(input_value)

if __name__=='__main__':
    app.run_server()
```

运行之后会的界面只有一个`dcc.Input`提供的输入框，但是这个输入框是输入后，是可以改变页面中的文字。那么这个是如何实现的呢？

我们的应用界面的输入和输出是通过`app.callback`装饰器进行声明。

在Dash中，应用的输入输出其实就是某个组件的属性(properties)。因此，`Output(component_id='my-div', component_property='children')`就可以解释为，将值输出到ID为`my-div`的HTML组件的`children`的参数中，而`[Input(component_id='my-id', component_property='value')]`则表明输入时来自于ID为`my-id`的`value`参数。

随着输入的值的改变，装饰器会调用函数`update_output_div`生成新值。这其实有点像Excel，当你写好一个函数后，修改原来值会产生新的值，这种编程方法叫做"Reactive Programming"，应该可以翻译为响应式编程吧.

让我们更进一步，看看使用`Slider`组件加上响应式编程后，图片是如何动起来. 数据和之前使用的一致，之前是展示了所有年份，不同洲的国家的GDP分布情况。而这里则可以使用滑动栏的方式，逐年查看。

```Python
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import pandas as pd

df = pd.read_csv(
    'https://raw.githubusercontent.com/plotly/'
    'datasets/master/gapminderDataFiveYear.csv')

app = dash.Dash()

app.layout = html.Div([
    dcc.Graph(id = 'graph-with-slider'),
    dcc.Slider(
        id = 'years-slider',
        min = df['year'].min(),
        max = df['year'].max(),
        value = df['year'].min(),
        step = None,
        marks = {str(year): str(year) for year in df['year'].unique()}
    )
])

@app.callback(
    dash.dependencies.Output(component_id = 'graph-with-slider', component_property = "figure"),
    [dash.dependencies.Input('years-slider', 'value')]
)
def update_figure(selected_year):
    filtered_df = df[df.year == selected_year]
    traces = []
    for i in filtered_df.continent.unique():
        df_by_continent = filtered_df[filtered_df['continent'] == i]
        traces.append(go.Scatter(
            x = df_by_continent['gdpPercap'],
            y = df_by_continent['lifeExp'],
            text = df_by_continent['country'],
            mode = 'markers',
            opacity = 0.7,
            marker = {
                'size': 15,
                'line': {'width':0.5, 'color':'white'}
            },
            name = i
        ))
    return {
        'data': traces,
        'layout': go.Layout(
            xaxis = {'type':'log', 'title':'GDP Per Capita'},
            yaxis = {'title':'Life Expectancy', 'range':[20,90]},
            margin = {'l':40, 'b':40, 't':10, 'r':10},
            legend = {'x':0, 'y':1},
            hovermode = 'closest'
        )
    }

if __name__ == '__main__':
    app.run_server()
```

首先是在布局中设置了两个占位组件，这两个占位组件一个用于提供年份用于筛选，一个用于则是展示输出。然后`update_figure`接受值返回对应的图形对象，最后展示到浏览器中。

Dash应用在启动的时候会加载数据，因此当用户访问应用的时候，数据已经在内存中，随后用户的交互操作就能得到及时的响应。当然`callback`函数不会修改原始数据，它仅仅是在内存中创建新的拷贝而已。

### 多个输入值

上一节只是单个输入单个输出，在Dash中，每个Output，都可以由多个Input。这一部分则是介绍通过加入更多调节组件多角度地展示数据。这里用到了五个调节组件，为2个`Dropdown`, 2个`RadioItems`和1个`Slider`。

```Python
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import pandas as pd

app = dash.Dash()

df = pd.read_csv(
    'https://gist.githubusercontent.com/chriddyp/'
    'cb5392c35661370d95f300086accea51/raw/'
    '8e0768211f6b747c0db42a9ce9a0937dafcbd8b2/'
    'indicators.csv')

available_indicators = df['Indicator Name'].unique()

app.layout = html.Div([
    html.Div([
        html.Div([
            dcc.Dropdown(
                id='xaxis-column',
                options=[{'label':i, 'value':i} for i in available_indicators],
                value = 'Fertility rate, total(births per woman)'
                ),
            dcc.RadioItems(
                id = 'xaxis-type',
                options = [{'label':i, 'value':i} for i in ['Liner','Log']],
                value = 'Liner',
                labelStype={'display':'inline-block'}
            )
        ],
        style = {'width':'48%', 'display':'inline-block'}),
        html.Div([
            dcc.Dropdown(
                id = 'yaxis-column',
                options = [{'label':i, 'value':i} for i in available_indicators],
                value = 'Life expectancy at birth, total(year)'
            ),
            dcc.RadioItems(
                id = 'yaxis-type',
                options = [{'label':i, 'value':i} for i in ['Liner','Log']],
                value = 'Liner',
                labelStyle={'display':'inline-block'}
            )
        ], style={'width':'48%','float':'right','display':'inline-block'})
    ]),
    dcc.Graph(id='indicator-graphic'),
    dcc.Slider(
        id='year-slider',
        min=df['Year'].min(),
        max=df['Year'].max(),
        value=df['Year'].max(),
        step=None,
        marks={str(year): str(year) for year in df['Year'].unique()}
    )
])

@app.callback(
    Output('indicator-graphic','figure'),
    [Input('xaxis-column','value'),
     Input('yaxis-column','value'),
     Input('xaxis-type','value'),
     Input('yaxis-type','value'),
     Input('year-slider','value')
    ]
)
def update_graph(xaxis_column_name, yaxis_column_name,
                 xaxis_type, yaxis_type,
                 year_value):
    dff = df[df['Year'] == year_value]
    return {
        'data':[go.Scatter(
            x=dff[dff['Indicator Name'] == xaxis_column_name]['Value'],
            y=dff[dff['Indicator Name'] == yaxis_column_name]['Value'],
            text=dff[dff['Indicator Name'] == yaxis_column_name]['Country Name'],
            mode = 'markers',
            marker = {
                'size': 15,
                'opacity': 0.5,
                'line':{'width':0.5, 'color':'white'}
            }
        )],
        'layout':go.Layout(
            xaxis={
                'title':xaxis_column_name,
                'type':'linear' if xaxis_type == 'Liner' else 'log'
            },
            yaxis={
                'title': yaxis_column_name,
                'type': 'linear' if yaxis_type == 'Liner' else 'log'
            },
            margin={'l':40, 'b':40,'t':10,'r':0},
            hovermode='closest'
        )
    }

if __name__ == '__main__':
    app.run_server()
```

和单个输入区别不大，就是输入多了，要写的代码多了，写代码的时候可能会写错而已。如果有多个输出的需求，只要定义多个`callback`函数即可。

### 第二部分小结

Dash应用使用装饰器`callback`进行响应式编程。回调函数根据`component_id`和`component_property`从不同组件中获取输入值，然后其所装饰的函数进行计算后，将值返回装饰器，最后将计算结果输出到指定组件中。