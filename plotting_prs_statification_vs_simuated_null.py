#!/usr/bin/env python

from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
import numpy as np
import rpy2
import base64

get_ipython().run_line_magic('load_ext', 'rpy2.ipython')

init_notebook_mode(connected=True)

# define input files
simulated_qvals_fh = "test_topGO_10K_raw_qvals.txt"
simulated_sample_dict_fh = "test_topGO_10K_sample_dict.txt"
simulated_gene_dict_fh = "test_topGO_10K_gene_table.txt"
observed_data_fh = '../first_EA_p1_data.txt'
observed_genes_fh = '../first_EA_p1_gene_list.txt'
observed_samples_fh = '../first_EA_p1_sample_list.txt'

## load data and reformat as necessary
# load results from simulation
null_qvals = pd.read_table(simulated_qvals_fh)
null_sample_dict = pd.read_table(simulated_sample_dict_fh)
null_gene_dict = pd.read_table(simulated_gene_dict_fh)

# load results from PRS stratification
# stats from topGO
observed_df = pd.read_table(observed_data_fh)
observed_df.columns = ['GO_term', 'Name', 'Annotated', 'Significant', 'Expected', 'p_val', 'q_val']

# gene list used in topGO
observed_genes = pd.read_table(observed_genes_fh)
observed_genes.columns = ['observed']
observed_genes = observed_genes.drop_duplicates()

# samples that have those genes
observed_samples = pd.read_table(observed_samples_fh)
observed_samples.columns = ['observed']

# -log 10 conversion
nlog_null_qvals = null_qvals.apply(lambda x: -np.log(x), axis = 0)
observed_df['nlog_q_val'] = -np.log(observed_df.q_val)
observed_df['nlog_p_val'] = -np.log(observed_df.p_val)

# build dict of GO term to name
go_dict = pd.Series(observed_df.Name.values, index = observed_df.GO_term).to_dict()

## new data frame for scatter plot
scatter_df = observed_df.loc[observed_df.q_val <= 0.05]

# get matching simulated values and reorder to match the observed values
null_scatter_df = nlog_null_qvals.loc[:, nlog_null_qvals.columns.isin(list(scatter_df.GO_term))]
null_scatter_df = null_scatter_df[scatter_df.GO_term]

## plot mean and SD of null vs PRS 
trace0 = go.Scatter(
    x = list(range(1,scatter_df.shape[0]+1)),
    y = null_scatter_df.mean(),
    marker = {'color': 'blue'},
    mode = 'markers',
    name = 'Simulated Null',
    error_y = dict(
        array = null_scatter_df.std()*2,
        visible = True),
    text = scatter_df.Name
)

trace1 = go.Scatter(
    x = list(range(1,observed_df.shape[0]+1)),
    y = scatter_df.nlog_q_val,
    marker = {'color': 'red'},
    mode = 'markers',
    name = 'PRS Stratification',
    text = scatter_df.GO_term
)

data = [trace0, trace1]

layout = dict(title = 'PRS stratification vs proband ascertainment bias',
              yaxis = dict(title = '-log Q-value (mean and 2 SDs)'),
              xaxis = dict(title = 'GO Term',ticks='',showticklabels=False)
)


fig = dict(data=data, layout=layout)
stick_plot_div = plot(fig, validate = False, include_plotlyjs=False, output_type='div')

## plot fold enrichment as a function of -log p-val
# color by q-value
trace0 = go.Scatter(
    x = observed_df.loc[observed_df.q_val > 0.05].nlog_p_val,
    y = observed_df.loc[observed_df.q_val > 0.05].Significant/observed_df.loc[observed_df.q_val > 0.05].Expected,
    marker = {'color': 'magenta'},
    mode = 'markers',
    name = 'q > 0.05',
    text = observed_df.loc[observed_df.q_val > 0.05].Name
)

trace1 = go.Scatter(
    x = observed_df.loc[observed_df.q_val < 0.05].nlog_p_val,
    y = observed_df.loc[observed_df.q_val < 0.05].Significant/observed_df.loc[observed_df.q_val < 0.05].Expected,
    marker = {'color': 'green'},
    mode = 'markers',
    name = 'q < 0.05',
    text = observed_df.loc[observed_df.q_val < 0.05].Name
)

data = [trace0, trace1]

layout = dict(title = 'Fold enrichment vs -log(p-value)',
            yaxis = dict(title = 'Fold enrichment (observed/expected)'),
            xaxis = dict(title = '-log(p-value)'),
            hovermode = 'closest'
)

fig = dict(data=data, layout=layout)

#iplot(fig, validate = False)
enrichment_by_pval_div = plot(fig, validate = False, include_plotlyjs=False, output_type='div')

## plot fold enrichment as a function of number of genes in GO term
# color by q-value
trace0 = go.Scatter(
    x = observed_df.loc[observed_df.q_val > 0.05].Annotated,
    y = observed_df.loc[observed_df.q_val > 0.05].Significant/observed_df.loc[observed_df.q_val > 0.05].Expected,
    marker = {'color': 'orange'},
    mode = 'markers',
    name = 'q > 0.05',
    text = observed_df.loc[observed_df.q_val > 0.05].Name
)

trace1 = go.Scatter(
    x = observed_df.loc[observed_df.q_val < 0.05].Annotated,
    y = observed_df.loc[observed_df.q_val < 0.05].Significant/observed_df.loc[observed_df.q_val < 0.05].Expected,
    marker = {'color': 'purple'},
    mode = 'markers',
    name = 'q < 0.05',
    text = observed_df.loc[observed_df.q_val < 0.05].Name
)

data = [trace0, trace1]

layout = dict(title = 'Fold enrichment vs GO term size',
            yaxis = dict(title = 'Fold enrichment (observed/expected)'),
            xaxis = dict(title = 'Number of Genes in GO term'),
            hovermode = 'closest'
)

fig = dict(data=data, layout=layout)

#iplot(fig, validate = False)
enrichment_by_size_div = plot(fig, validate = False, include_plotlyjs=False, output_type='div')

## plotting histogram of null q-values
test_go_term = "GO:0086010" # user input goes here

x = nlog_null_qvals[test_go_term]

trace0 = go.Histogram(
    x = x,
    nbinsx = 100,
    text = nlog_null_qvals.index
)

data = [trace0]

layout = {'title': go_dict[test_go_term], 
    'yaxis': dict(title = 'Frequency'),
    'xaxis': dict(title = '-log(Simulated q-value) vs. -log(observed q-value)'),
    'hovermode': 'closest',
    'bargap': 0.1, 
    'shapes':[{'type': 'line',
        'x0':observed_df.loc[observed_df.GO_term == test_go_term].nlog_q_val.values[0], 
        'y0':0,
        'x1':observed_df.loc[observed_df.GO_term == test_go_term].nlog_q_val.values[0],
        'y1':25,
        'line':{'color':'red','width': 3}
    }]
}

fig = dict(data=data, layout=layout)

#iplot(fig, validate = False)
hist_div = plot(fig, validate = False, include_plotlyjs=False, output_type='div')

## Explore simulations with test statistics more extreme than observed statistic
## in prep for plottin venn diagrams below

# go term of interest
test_go_term = "GO:0086010" # user input goes here

# observed value
test_stat = observed_df.loc[observed_df.GO_term == test_go_term].nlog_q_val.values[0]

# select simulation data which are more extreme
outlier = nlog_null_qvals.loc[nlog_null_qvals[test_go_term] > test_stat].index.values

# lookup in gene list dictionary
outlier_samples = null_sample_dict[outlier]

# lookup in sample dictionary
outlier_genes = null_gene_dict.loc[null_gene_dict.index.isin(outlier)]
outlier_genes = outlier_genes.loc[:, (outlier_genes != 0).any(axis = 0)]
#outlier_genes = outlier_genes.transpose()

# add on the observed data
# samples
outlier_samples['observed'] = observed_samples.observed.values

# genes
outlier_genes = outlier_genes.transpose()
outlier_genes['observed'] = 0

for i in outlier_genes.index.values:
    if i in observed_genes.observed.values:
        outlier_genes.at[i,'observed'] = 1
    else:
        continue

# write to file for potting
outlier_samples.to_csv("test_outlier_samples.txt", sep = "\t", index = False)
outlier_genes.to_csv("test_outlier_genes.txt", sep = "\t", index = True)

get_ipython().run_cell_magic('R', '', '## plot venn diagaram of overlap among samples in groups\n\n# suppress venn.diagram logging\nfutile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")\n\n# sample overlap\nlibrary(VennDiagram)\noutlier_samples = read.table("test_outlier_samples.txt", header = TRUE)\nhead(outlier_samples)\n\nx = list()\nfor(i in 1:ncol(outlier_samples)){\n    x[[i]] = outlier_samples[,i]\n}\n\nvenn.diagram(x, \n    category.names = colnames(outlier_samples),\n    filename = "test_samples.jpeg",\n    output = TRUE,\n    cex = 0.75,\n    cat.cex = 0.5\n)')

get_ipython().run_cell_magic('R', '', '## plot venn diagram of overlap among genes in groups\n\n# suppress venn.diagram logging\nfutile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")\n\n# gene overlap\nlibrary(VennDiagram)\noutlier_genes = read.table("test_outlier_genes.txt", header = TRUE)\n\nx = list()\nfor(i in 1:ncol(outlier_genes)){\n    y = c()\n    for(j in 1:nrow(outlier_genes)){\n        if(outlier_genes[j,i] == 0){\n            next\n        }\n        else {\n            y = c(y, row.names(outlier_genes)[j])\n        }\n    }\n    x[[i]] = y    \n}\n\nvenn.diagram(x,\n    category.names = colnames(outlier_samples),\n    filename = "test_genes.png",\n    cex = 0.75,\n    cat.cex = 0.5\n)')

## package the results
# base64 encode images
samples_img = base64.b64encode(open("test_samples.png", "rb").read()).decode('utf-8')
genes_img = base64.b64encode(open("test_genes.png", "rb").read()).decode('utf-8')

# build an html string that contains the plots above
html_string = '''
<html>
    <!-- source plotly js -->
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    
    <head>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
        <style>body{ margin:0 100; background:whitesmoke; }</style>
    </head>
    
    <body>
        <h1>First Quartile of PRS for EA in Probands</h1>

        <!-- *** Section 1 *** --->
        <h2>Simulated vs observed enrichment of GO terms</h2>
        '''+stick_plot_div+'''
        <p>Red dots outside of blue bars are GO terms significant beyond proband ascertainment bias. Hover over \
        lines to get GO term ID, GO term name, observed -log q-value, and simulated mean and standard deviation.</p>
        
        <!-- *** Section 2 *** --->
        <h2>Fold enrichment as a function of p-value</h2>
        '''+enrichment_by_pval_div+'''
        <p>Fold enrichment is calculated as the number of observed genes divided by the number of expected. \
        Hover over points to get GO term ID.</p>

        <!-- *** Section 3 *** --->
        <h2>Fold enrichment as a function total genes in GO term</h2>
        '''+enrichment_by_size_div+'''
        <p>Fold enrichment is calculated as the number of observed genes divided by the number of expected. \
        Hover over points to get GO term ID.</p>

        <!-- *** Section 4 *** --->
        <h2>Histogram of GO terms that are more significant than expected by proband ascertainment bias</h2>
        '''+hist_div+'''
        <p>Vertical red line indicates value of observed test statistic. Hover over bars to get simulation IDs from \
        Monte Carlo used to generate null hypothesis.</p>

        <!-- *** Section 5 *** --->
        <h2>Venn diagram of sample overlap</h2>
        <img src = "data;image/png;base64,'''+samples_img+'''", alt = "genes venn", width = "400", height = "400" >
        <p>Sample overlap between observed (red line in histogram) and simulations greater than this value.</p>

        <!-- *** Section 6 *** --->
        <h2>Venn diagram of gene overlap</h2>
        <img src = "data;image/png;base64,'''+genes_img+'''", alt = "genes venn", width = "400", height = "400" >
        <p>Sample overlap between genes used in observed (red line in histogram) and genes used in simulations \
        greater than this value.</p>

    </body>
</html>'''

# write to file
f = open('./test_report.html','w')
f.write(html_string)
f.close()
