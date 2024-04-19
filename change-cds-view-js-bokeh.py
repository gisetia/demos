# %%
from bokeh.io.output import output_file
from bokeh.plotting import figure, show, output_notebook
from bokeh.models import MultiSelect, BooleanFilter, Slider, CheckboxGroup, CustomJS, ColumnDataSource, CDSView
from bokeh.models.filters import CustomJSFilter
from bokeh.layouts import row
from bokeh.transform import factor_cmap
from bokeh.palettes import Category10_10
from bokeh.models import AutocompleteInput

from bokeh.sampledata import iris

source = ColumnDataSource(data=iris.flowers)
species = iris.flowers.species.unique().tolist()

checkboxes = CheckboxGroup(labels=species, active=list(range(len(species))), 
name='Species')

fig = figure(x_range=(4, 8), y_range=(1.5, 4.5))

bool_filter = BooleanFilter([s in species for s in source.data['species']])
view = CDSView(source=source, filters=[bool_filter])

fig.scatter("sepal_length", "sepal_width", source=source,
            color=factor_cmap("species", Category10_10, species), view=view)

checkboxes.js_on_change('active', CustomJS(
    args=dict(b=bool_filter, s=source, labels=species),
    code="""
        let active_ind = this.active;
        let active_elem = active_ind.map(index => labels[index]);
        let column_name = this.name.toLowerCase()
        b.booleans = s.data[column_name].map(element => active_elem.includes(element));
        s.change.emit();
        console.log(this.name + ' ' + active_elem + ' ');

        """))

output_file('tmp-files/filter-data.html')
show(row(checkboxes, fig))

# %%
