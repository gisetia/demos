#%%
import numpy as np
import pandas as pd

from bokeh.plotting import figure, output_file, show
from bokeh.layouts import gridplot, column, row
from bokeh.models import (ColumnDataSource, HoverTool, Slope, Whisker, CDSView,
                          LinearColorMapper, LogColorMapper, CheckboxGroup)

# corr = pd.read_parquet('processed_data/corr.pq')
sem_sim = pd.read_parquet('processed_data/sem-sim.pq')

# %%
