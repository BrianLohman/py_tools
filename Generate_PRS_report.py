#!/usr/bin/env python
# coding: utf-8

# In[1]:


from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
import numpy as np
import json
import re


# In[2]:


# define input and output files
simulated_qvals_fh = "test_topGO_10K_raw_qvals.txt"
simulated_sample_dict_fh = "test_topGO_10K_sample_dict.txt"
simulated_gene_dict_fh = "test_topGO_10K_gene_table.txt"
observed_data_fh = '../first_EA_p1_data.txt'
observed_genes_fh = '../first_EA_p1_gene_list.txt'
observed_samples_fh = '../first_EA_p1_sample_list.txt'

report = "test_report.html"


# In[3]:


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

#build a dict of name to GO term
name_dict = pd.Series(scatter_df.GO_term.values, index = scatter_df.Name).to_dict()


# In[4]:


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


# In[5]:


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

enrichment_by_pval_div = plot(fig, validate = False, include_plotlyjs=False, output_type='div')


# In[6]:


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

enrichment_by_size_div = plot(fig, validate = False, include_plotlyjs=False, output_type='div')


# In[7]:


# master dict to hold all data relvant to each GO term
master_GO_dict = {}

# list of GO names for dropdown menus later
names = []

all_go_terms = scatter_df.GO_term
#all_go_terms = ["GO:0086010"]

for test_go_term in all_go_terms: 
    # start dictionary with basic info
    GO = {}
    GO['go_id'] = test_go_term
    GO['name'] = go_dict[test_go_term] # lookup from prior dict 

    # save list of names for dropdown menus later
    names.append(go_dict[test_go_term])
    
    ## plotting histogram of null q-values
    GO['hist_vals'] = list(nlog_null_qvals[test_go_term].values)
    GO['hist_line'] = observed_df.loc[observed_df.GO_term == test_go_term].nlog_q_val.values[0]
    
    ## Explore simulations with test statistics more extreme than observed statistic
    # observed value
    test_stat = observed_df.loc[observed_df.GO_term == test_go_term].nlog_q_val.values[0]

    # select simulation data which are more extreme
    outlier = nlog_null_qvals.loc[nlog_null_qvals[test_go_term] > test_stat].index.values

    # lookup in gene list dictionary
    outlier_samples = null_sample_dict[outlier]

    # lookup in sample dictionary
    outlier_genes = null_gene_dict.loc[null_gene_dict.index.isin(outlier)]
    outlier_genes = outlier_genes.loc[:, (outlier_genes != 0).any(axis = 0)]

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

    outlier_genes['Total_Observations'] = outlier_genes.sum(axis = 1)
    outlier_genes = outlier_genes[outlier_genes.columns[::-1]]
    #print(outlier_genes.head())

    # convert to list of lists (where each list is a row) for jquery
    gene_header = ['gene']
    gene_header.extend(outlier_genes.columns)

    #print(header)
    gene_header = [{"title":h} for h in gene_header]

    # save to dict
    GO['gene_header'] = gene_header
    
    gene_table = []
    for i in range(len(outlier_genes.index)):
        gene = outlier_genes.index.values[i]
        tmp = [gene]
        tmp.extend(outlier_genes.iloc[i,].values)
        tmp = [str(x) for x in tmp]
        gene_table.append(tmp)

    # add to dictionary
    GO['gene_table'] = gene_table

    ## outlier samples
    #print(outlier_samples.head())

    df_arr = outlier_samples.values
    groups = outlier_samples.columns.values
    colnames = ['sample_id']
    colnames.extend(groups)

    empty_df = pd.DataFrame(columns=colnames)

    for i, sample in enumerate(set(df_arr.flatten())):
        tf_array = sample == df_arr
        a = (np.array(np.where(tf_array == True)))
        ii = a[-1]
        vals = np.zeros(len(groups))
        vals[ii] = 1
        vals = list(vals)
        col = [sample]
        vals.insert(0, sample)
        empty_df = empty_df.append(dict(zip(colnames, vals)), ignore_index = True)

    sample_table = empty_df
    sample_table['Total_Observations'] = sample_table.sum(axis = 1)

    sample_table = sample_table.set_index(sample_table.sample_id)
    sample_table.drop('sample_id', inplace = True, axis = 1)
    sample_table = sample_table.astype(int)
    sample_table = sample_table[sample_table.columns[::-1]]
    #print(sample_table.head())

    # convert to list of lists (where each list is a row) for jquery
    sample_header = ['sample']
    sample_header.extend(outlier_genes.columns)

    #print(header)
    sample_header = [{"title":h} for h in sample_header]

    # save to dict
    GO['sample_header'] = sample_header
    
    sample_lists = []
    for i in range(len(outlier_samples.index)):
        sample = sample_table.index.values[i]
        tmp = [sample]
        tmp.extend(sample_table.iloc[i,].values)
        tmp = [str(x) for x in tmp]
        sample_lists.append(tmp)

    # add to dictionary
    GO['sample_table'] = sample_lists
    
    # add all data from each GO term to dict, keyed by GO term ID
    master_GO_dict[test_go_term] = GO


# In[62]:


TEMPLATE = '''
<!DOCTYPE html>



<html>
    <head>
        <meta charset="utf-8"/>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <script type="text/javascript" src="https://code.jquery.com/jquery-3.3.1.js"></script>
        <script type="text/javascript" src="https://cdn.datatables.net/v/bs4/dt-1.10.18/sl-1.3.0/datatables.min.js"></script>
        <script src="https://cdn.datatables.net/1.10.16/js/dataTables.bootstrap4.min.js" type='text/javascript'></script>
        <link href="https://cdn.datatables.net/v/bs4/dt-1.10.18/sl-1.3.0/datatables.min.css" rel="stylesheet" type="text/css"/>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">

        <style>
            body {margin: 1% 2.5%; background:white }
            p { margin: 0 2.5% }
            h2 { margin: 1% 0}
            .dataTables_filter {
                text-align: left !important;
            }
        </style>
    </head>
    
    <body>
        <h1 align = 'center'>First Quartile of PRS for EA in Probands</h1>

        <!- - *** Section 1 *** - ->
        <h2>Simulated vs observed enrichment of GO terms</h2>
        [STICK_PLOT_DIV]
        <p>Red dots outside of blue bars are GO terms significant beyond proband ascertainment bias. Hover over \
        lines to get GO term ID, GO term name, observed -log q-value, and simulated mean and standard deviation.</p>
        
        <!- - *** Section 2 *** - ->
        <h2>Fold enrichment as a function of p-value</h2>
        [PVAL_DIV]
        <p>Fold enrichment is calculated as the number of observed genes divided by the number of expected. \
        Hover over points to get GO term ID.</p>

        <!- - *** Section 3 *** - ->
        <h2>Fold enrichment as a function total genes in GO term</h2>
        [SIZE_DIV]
        <p>Fold enrichment is calculated as the number of observed genes divided by the number of expected. \
        Hover over points to get GO term ID.</p>

        <!- - *** Section 4 *** - ->
        <h2> Analysis of Specific GO terms (populates histogram and tables below) </h2>
        <div>
            <label>Select GO Term</label>
            <select id="go_select">
                [TERM_SELECT]
            </select>
        </div>
        
        <h2>Histogram of GO terms that are more significant than expected by proband ascertainment bias</h2>
        <div id = 'hist_div_id'> </div>
        
        <p>Vertical red line indicates value of observed test statistic. Hover over bars to get simulation IDs from \
        Monte Carlo used to generate null hypothesis.</p>

        <!- - *** Section 5 *** - ->
        <h2>Genes in simulated and observed groups</h2>
        <p>Cells indicate presence (1) or absence (0) of a gene in a group. Click colums to sort. \
        Genes present in all groups may represent a core set of genes driving GO enrichment \
        (Sort based on Total Observations). Control + click to sort based on second category. </p>
        <table id="gene_table" class="table table-hover pb-3 display nowrap" width="90%"> </table>

        <h2>Samples in simulated and observed groups</h2>
        <p>Cells indicate presence (1) or absence (0) of a sample in a group. Click colums to sort. \
         Control + click to sort based on second category.</p>
        <table id="sample_table" class="table table-hover pb-3 display nowrap" width="90%"></table>

    </body>

<script>
const data = [DATA]
var gene_table
var sample_table
var go = Object.keys(data)[0]

const build_hist = (go) => {
    hist_vals = data[go]['hist_vals']
    hist_line = data[go]['hist_line']
    hist_title = data[go]['name']
    
    var trace = {
        x: hist_vals,
        type: 'histogram'
    }
    var hist_data = [trace]
    
    hist_layout = {'title': hist_title, 
        'yaxis': {title: 'Frequency'},
        'xaxis': {title: '-log(Simulated q-value) vs. -log(observed q-value)'},
        'hovermode': 'closest',
        'bargap': 0.1, 
        'shapes':[{'type': 'line',
            'x0':hist_line, 
            'y0':0,
            'x1':hist_line,
            'y1':1,
            'yref': "paper",
            'line':{'color':'red','width': 3}
        }]
    }
    
    let hist_plot = document.getElementById("hist_div_id")
    Plotly.react(hist_plot, hist_data, hist_layout)

}

const build_table = (go) => {
    if ( $.fn.DataTable.isDataTable('#gene_table') ) {
      $('#gene_table').DataTable().destroy()
    }
    $('#gene_table tbody').empty()
    if ( $.fn.DataTable.isDataTable('#sample_table') ) {
      $('#sample_table').DataTable().destroy()  
    }
    $('#sample_table tbody').empty()
    
    gene_table = $("#gene_table").DataTable({
        data: data[go]['gene_table'],
        columns: data[go]['gene_header'], 
        scrollY: '600px',
        scrollX: true,
        scrollCollapse: true,
        paging: false,
        pagingType: "simple",
        info: true,
        dom: 'flrtip'
    })

    sample_table = $("#sample_table").DataTable({
        data: data[go]['sample_table'],
        columns: data[go]['sample_header'], 
        scrollY: '600px',
        scrollX: true,
        scrollCollapse: true,
        paging: false,
        pagingType: "simple",
        info: true,
        dom: 'flrtip'
    })
}

$('#go_select').on("change", () => { 
    go = $('#go_select :selected').val()
    build_table(go)
    build_hist(go)
})

$(document).ready(function() {
    build_table(go)
    build_hist(go)
})

</script>
</html>
'''


# In[63]:


# add in the data
html = TEMPLATE.replace("[DATA]", json.dumps(master_GO_dict))

# plots that only depend on the focal PRS, quartile, and proband/sibling status
html = html.replace("[STICK_PLOT_DIV]", str(stick_plot_div))
html = html.replace("[PVAL_DIV]", str(enrichment_by_pval_div))
html = html.replace("[SIZE_DIV]", str(enrichment_by_size_div))

select_term = []
for key in name_dict.keys():
    select_term.append("<option value = "+name_dict[key]+">"+key+"</option>")
select_term = " ".join(select_term)

html = html.replace("[TERM_SELECT]", select_term)

f = open(report, "w")
f.write(html)
f.close()


# In[ ]:




