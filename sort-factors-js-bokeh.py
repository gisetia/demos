# %%
import pandas as pd
from bokeh.io import output_file, show
from bokeh.plotting import figure
from bokeh.models import MultiSelect, CustomJS, ColumnDataSource, RadioGroup
from bokeh.layouts import row, column


output_file("test.html")

fruits = ['Apples', 'Pears', 'Nectarines', 'Plums', 'Grapes',
          'Strawberries']
counts = [5, 3, 4, 2, 4, 6]
counts2 = [8, 7, 5, 3, 2, 1]

src = ColumnDataSource({'fruits': fruits, 'counts': counts,
                        'counts2': counts2})


p = figure(x_range=fruits, plot_height=250, title="Fruit Counts",
           toolbar_location=None, tools="")

p.vbar(x='fruits', top='counts', width=0.9, source=src)

p.xgrid.grid_line_color = None
p.y_range.start = 0

# Select column to sort factors
select = MultiSelect(title="Select:", value=[''],
                     options=['counts', 'counts2'])


callback_code = """
var data = src.data;

var to_sort = select.value;


/*const a = data.fruits.map((v, i) => 
({fruit: v, sort: data[to_sort][i]}))
    .sort((a, b) => a.sort - b.sort).map((v) => v.fruit);*/

const a = data['fruits'].map((v, i) => 
({fruit: v, sort: data[to_sort][i]}))
    .sort((a, b) => a.sort - b.sort).map((v) => v['fruit']);

p.x_range.factors = a

console.log( '  ' + to_sort)
"""

callback = CustomJS(
    args=dict(src=src,
              select=select, p=p),
    code=callback_code)

select.js_on_change('value', callback)

# labels = [['Counts', 'Counts 2'],
# ('counts', 'counts2')]

radio_group = RadioGroup(labels=labels[0], active=0)

# callback_code_radio = '''
# var data = src.data;
# var to_sort = labels[1][radio.active]

# const a = data.fruits.map((v, i) => ({fruit: v, sort: data[to_sort][i]}))
#     .sort((a, b) => a.sort - b.sort).map((v) => v.fruit);


# p.x_range.factors = a


# console.log(to_sort)
# '''

# callback_radio = CustomJS(
#     args=dict(src=src,
#               radio=radio_group, p=p, labels=labels),
#     code=callback_code_radio)

# radio_group.js_on_click(callback_radio)

show(row(select, radio_group, p))

# %%
