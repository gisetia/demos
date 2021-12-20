# %%
from bokeh.io.output import output_file
from bokeh.plotting import figure, show, output_notebook
from bokeh.models import Slider, CheckboxGroup, CustomJS, ColumnDataSource, CDSView
from bokeh.models.filters import CustomJSFilter
from bokeh.layouts import row
from bokeh.transform import factor_cmap
from bokeh.palettes import Category10_10
from bokeh.models import AutocompleteInput

from bokeh.sampledata import iris

source = ColumnDataSource(data=iris.flowers)
species = iris.flowers.species.unique().tolist()

checkboxes = CheckboxGroup(labels=species, active=list(range(len(species))))
autocomp = AutocompleteInput(value='', completions=species, min_characters=1)

fig = figure(x_range=(4, 8), y_range=(1.5, 4.5))

filter = CustomJSFilter(code="""
//let selected = checkboxes.active.map(i=>checkboxes.labels[i]);

var selected = autocomp.value;
var indices = [];
var column = source.data.species;
for(var i=0; i<column.length; i++){
    if(selected.includes(column[i])){
        indices.push(i);
    }
}
return indices;
cosole.log('aaa' + selected)
""", args=dict(checkboxes=checkboxes, autocomp=autocomp))

checkboxes.js_on_change("active", CustomJS(code="source.change.emit();",
                                           args=dict(source=source)))
autocomp.js_on_change("value", CustomJS(code="source.change.emit();",
                                        args=dict(source=source)))

fig.scatter("sepal_length", "sepal_width",
            color=factor_cmap("species", Category10_10, species),
            source=source, view=CDSView(source=source, filters=[filter]))
output_file('tmp-files/filter-data.html')
show(row(autocomp, checkboxes, fig))

# %%
