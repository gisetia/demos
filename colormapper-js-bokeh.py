import pandas as pd
from bokeh.models import ColumnDataSource, ColorBar, Select, CustomJS
from bokeh.plotting import figure
from bokeh.layouts import gridplot
from bokeh.palettes import Spectral5, Viridis5
from bokeh.transform import linear_cmap
from bokeh.embed import components
from jinja2 import Environment, FileSystemLoader


df = pd.DataFrame({"a": [2, 6, 5, 3, 7, 8, 1, 9, 2, 4],
                   "b": [3, 5, 7, 1, 0, 6, 5, 4, 2, 9],
                   "c": [11, 12, 13, 14, 11, 13, 15, 14, 15, 12],
                   "d": [21, 23, 24, 25, 21, 22, 23, 24, 25, 22]})
source = ColumnDataSource(df)

mapper = linear_cmap(field_name="c", palette=Spectral5,
                     low=min(df["c"]), high=max(df["c"]))

fig = figure(plot_width=400, plot_height=400)

cir = fig.circle(x="a", y="b", size=12,
                 source=source, line_color=mapper, color=mapper)
color_bar = ColorBar(color_mapper=mapper["transform"], width=8,
                     location=(0, 0))
fig.add_layout(color_bar, "right")

codec = """
            var low = Math.min.apply(Math,source.data[cb_obj.value]);
            var high = Math.max.apply(Math,source.data[cb_obj.value]);
            var color_mapper = new Bokeh.LinearColorMapper({palette:pal, low:low, high:high});
            cir.glyph.fill_color = {field: cb_obj.value, transform: color_mapper};
            cir.glyph.line_color = {field: cb_obj.value, transform: color_mapper};
            //color_bar.color_mapper = color_mapper;
            color_bar.color_mapper.low = low;
            color_bar.color_mapper.high = high;
            color_bar.color_mapper.palette = pal;
            source.change.emit();
        """
cb_cselect_c = CustomJS(args=dict(cir=cir, source=source, color_bar=color_bar, pal=Viridis5),
                        code=codec)

c_select = Select(title="Select variable for color: ", value="None",
                  options=["c", "d"], callback=cb_cselect_c)

layout = gridplot([[fig], [c_select]])

env = Environment(loader=FileSystemLoader("."))
template = env.get_template("template.html")

script, div = components(layout)

with open("output.html", "w") as f:
    print(template.render(script=script, div=div), file=f)
