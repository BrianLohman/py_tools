#!/usr/bin/env python
# make a table of the number of impacts (low,med,high) for each individual in the variant table
# make a plot of the distribution of med and high impact variants
import sys

import numpy as np
import pandas
import plotly.graph_objs as go
import plotly.plotly as py
from plotly.offline import download_plotlyjs
from plotly.offline import init_notebook_mode
from plotly.offline import plot

init_notebook_mode(connected=True)

# define I/O
variants = pandas.read_table(sys.argv[1], low_memory=False)
out = open("individuals_by_impact_" + sys.argv[1], "w")

# count of low, med, and high impact mutations for each individual
# make a list of all the sample names
samples = list(variants)[8:]

# transform to individuals as rows, filtering any 0 and -1 counts
result = []
for i in samples:
    x = []
    x.append(i)
    tmp = variants[["impact", i]]
    x.append(tmp.loc[tmp[i] >= 1].impact.value_counts())
    result.append(x)

# reorganize to final output
data = dict()
for i in result:
    sample, series = i
    data[sample] = dict()
    for j in ["LOW", "MED", "HIGH"]:
        try:
            data[sample][j] = series[j]
        except KeyError:
            pass

# save observations of low and high, send to file
med = []
high = []
sample_id = []

for sample, values in data.iteritems():
    sample_id.append(sample)
    med.append(values["MED"])
    high_val = 0
    try:
        high_val = values["HIGH"]
    except KeyError:
        pass
    high.append(high_val)
    out.write(
        "\t".join(
            [sample, str(values["LOW"]), str(values["MED"]), str(high_val) + "\n"]
        )
    )

# make the plot
trace0 = go.Box(
    y=np.asarray(med), name="MED Impact", text=sample_id, hoverinfo="y+text"
)
trace1 = go.Box(
    y=np.asarray(high), name="HIGH Impact", text=sample_id, hoverinfo="y+text"
)
plot_data = [trace0, trace1]

layout = go.Layout(
    title="Variants by predicted impact",
    yaxis=dict(title="Number of Variants", zeroline=False),
)

fig = go.Figure(data=plot_data, layout=layout)

plot(fig, filename=sys.argv[1] + "_variant_plot.html", validate=False, auto_open=False)
