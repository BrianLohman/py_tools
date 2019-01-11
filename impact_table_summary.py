# make a table of the number of impacts (low,med,high) for each individual in the variant table
# make a plot of the distribution of med and high impact variants

from plotly.offline import download_plotlyjs, init_notebook_mode, plot
import pandas
import sys
import plotly.plotly as py
import plotly.graph_objs as go
import numpy as np
init_notebook_mode(connected=True)

# define I/O
variants = pandas.read_table(sys.argv[1], low_memory = False)
out = open("individuals_by_impact"+sys.argv[1], 'w')

# count of low, med, and high impact mutations for each individual
# make a list of all the sample names
samples = list(variants)[8:]

# transform to individuals as rows, filtering any 0 and -1 counts
result = []
for i in samples:
	x = []
	x.append(i)
	tmp = variants[['impact', i]]
	x.append(tmp.loc[tmp[i] >= 1].impact.value_counts())
	result.append(x)

# reorganize to final output
data = dict()
for i in result:
	sample, series = i
	data[sample] = dict()
	for j in ['LOW', 'MED', 'HIGH']:
		try:
			data[sample][j] = series[j]
		except KeyError:
			pass

# save observations of low and high, send to file
med = []
high = []
for sample, values in data.iteritems():
	med.append(values['MED'])
	high.append(values['HIGH'])
	out.write('\t'.join([sample, str(values['LOW']), str(values['MED']), str(values['HIGH'])]))

# make the plot
trace0 = go.Box(
	y = np.asarray(med)
)
trace1 = go.Box(
	y = np.asarray(high)
)
plot_data = [trace0, trace1]

layout = go.Layout(
	title = "Variants by predicted impact",
	yaxis=dict(title = 'Number of Variants', zeroline = False)
)

fig = go.Figure(data = plot_data, layout = layout)

plot(fig, filename = sys.argv[1]+"variant_plot.html", validate = False)
