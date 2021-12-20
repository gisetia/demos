# %%
import random
from bokeh.models import CustomJS, Div, ColumnDataSource
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row

x = [random.random() for x in range(500)]
y = [random.random() for x in range(500)]
txt = [str(random.randint(0, 100)) for x in range(500)]

source = ColumnDataSource(data=dict(x=x, y=y, txt=txt))
p = figure(tools='lasso_select,reset')
p.circle('x', 'y', source=source)

txt_lasso = Div(text='', margin=(10, 0, 0, 30), width=200)

code_lasso = '''
var inds = cb_obj.indices;
var txt = '';
for (var i = 0; i < inds.length; i++) {
            txt = txt + source.data['txt'][inds[i]] + '<br/>'
        }
txt_lasso.text = txt;
'''

lasso_callback = CustomJS(args=dict(source=source, txt_lasso=txt_lasso),
                          code=code_lasso)
source.selected.js_on_change('indices', lasso_callback)

output_file('txt_lasso.html')
layout = row(p, txt_lasso)
show(layout)


# %%
