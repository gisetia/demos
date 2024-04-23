#%%
import pandas as pd
from bokeh.plotting import figure, output_file, show
from bokeh.models import (ColumnDataSource, Range1d, HoverTool, 
                          LinearColorMapper, RangeTool)
from bokeh.palettes import Category10, Viridis256, Set1, Turbo256
from bokeh.layouts import row, column


umap_data = pd.read_csv('depmap_umap_data.csv')
clusters = pd.read_csv('depmap_umap_clusters.csv')

# umap_data.plot.scatter('x_coord', 'y_coord', s=2)

# %%

source = ColumnDataSource(umap_data.query('cluster_id >= 0'))
clust_source = ColumnDataSource(clusters.query('cluster_id != -1'))

zoomed = figure(width=500, height=500,
                x_range=Range1d(-15, -15,  bounds=(umap_data['x_coord'].min(),
                                                   umap_data['x_coord'].max())),
                y_range=Range1d(-15, -15, bounds=(umap_data['y_coord'].min(),
                                                  umap_data['y_coord'].max())),
                toolbar_location='below', tools='wheel_zoom, pan',
                margin=(10, 0, 0, 12))
zoomed.grid.visible = False
zoomed.xaxis.visible = False
zoomed.yaxis.visible = False
zoomed.toolbar.logo = None

mapper = LinearColorMapper(palette=list(reversed(Viridis256)),
                           low=-.1, high=1.01)
mapper_labels = LinearColorMapper(palette=Turbo256)

# zoomed.scatter(x='x_centroid', y='y_centroid', name='clusters',
#                source=clust_source, radius='plot_size',
#                alpha=0.5, fill_color='white',
#                line_color={'field': 'color', 'transform': mapper},
#                line_width=15, marker='circle', nonselection_line_color=None,
#                nonselection_alpha=0.5)
zoomed.scatter(x='x_coord', y='y_coord', source=source,
               name='genes',
            #    fill_color={'field': 'color', 'transform': mapper},
               line_color=None, size=25, alpha=0.5,
               selection_line_color='#e84393', selection_line_width=8,
               selection_alpha=1, nonselection_alpha=0.7, marker='circle')


zoomed.add_tools(HoverTool(tooltips=[('Gene', '@gene'),
                                     ('95% -log10(p)', '@perc'),
                                     ('Neg association ratio', '@sign'),
                                     ('In reactome term', '@in_reactome')],
                           names=['genes']))

plt = figure(width=800, height=800,
            #  title=(f'UMAP p_thr={p_thr}, eps={eps}, '
            #         f'min_dist_clust={min_dist_clust}, '
            #         f'min_dist_plot={min_dist_plot}, n_neighbors_clust='
            #         f'{n_neighbors_clust}, n_neighbors_plot='
            #         f'{n_neighbors_plot}'),
             toolbar_location='left'
             )
plt.xgrid.visible = False
plt.ygrid.visible = False
plt.xaxis.visible = False
plt.yaxis.visible = False
plt.toolbar.logo = None

range_tool = RangeTool(x_range=zoomed.x_range,
                       y_range=zoomed.y_range)
range_tool.overlay.fill_color = None
range_tool.overlay.line_color = None
plt.add_tools(range_tool)


# menus = column(row(gene_highlight, spinner), color_by)
# layout = row(plt, column(menus, zoomed), txt_cluster)

layout = row(plt, zoomed)

output_file(f'depmap-umap.html')
show(layout)
# %%
