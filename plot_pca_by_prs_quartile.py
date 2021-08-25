#!/usr/bin/env python

from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import pandas as pd
import matplotlib.pyplot as plt
import plotly.plotly as py
import plotly.graph_objs as go

'file:/mnt/c/Users/Lohman/Documents/Autism/temp-plot.html'
init_notebook_mode(connected=True)

df = pd.read_table("15_Jan_19_Simons_master_ancestry_corrected_PRS.txt")

# function to assign bins
def which_quartile(i):
    if i <= summary[4]:
        return "1"
    if i >= summary[6]:
        return "4"
    return "2_3"

moms = df.loc[df.family_member == "mo"]
summary = moms.NewBMI_ancestry_resid.describe()
moms['mom_bmi_quartile'] = moms['NewBMI_ancestry_resid'].apply(which_quartile)

#print(moms.head())
bmi_dict = {}
for fam in moms.family:
    bmi_dict[int(fam)] = moms.loc[moms.family == fam].mom_bmi_quartile.values[0]

probands = df.loc[df.family_member == "p1"]

x = []
for fam in probands.family:
#print(fam)
#print(probands.loc[probands.family == fam])
x.append(bmi_dict[fam])

probands['mom_bmi_quartile'] = x

summary = probands['EA_ancestry_resid'].describe()
probands['EA_quartile'] = probands['EA_ancestry_resid'].apply(which_quartile)

summary = probands['SCZ_ancestry_resid'].describe()
probands['SCZ_quartile'] = probands['SCZ_ancestry_resid'].apply(which_quartile)

trace0 = go.Scatter(
    x = probands.loc[probands.EA_quartile == '1'].plink_pc1,
    y = probands.loc[probands.EA_quartile == '1'].plink_pc2,
    marker = {'color': "red"},
    mode = 'markers',
    name = 'First Quartile',
    opacity = 0.75,
    text = probands.IID
)

trace1 = go.Scatter(
    x = probands.loc[probands.EA_quartile == '4'].plink_pc1,
    y = probands.loc[probands.EA_quartile == '4'].plink_pc2,
    marker = {'color': "blue"},
    mode = 'markers',
    name = 'Fourth Quartile',
    opacity = 0.75,
    text = probands.IID
)

trace2 = go.Scatter(
    x = probands.loc[probands.EA_quartile == '2_3'].plink_pc1,
    y = probands.loc[probands.EA_quartile == '2_3'].plink_pc2,
    marker = {'color': "grey"},
    mode = 'markers',
    name = 'Middle Quartiles',
    opacity = 0.5,
    text = probands.IID
)

data = [trace0, trace1, trace2]

layout = dict(title = 'Residual Ancestry PCA colored by EA PRS stratification',
    yaxis = dict(title = 'PC 2'),
    xaxis = dict(title = 'PC 1')
                                                 )

fig = dict(data=data, layout=layout)

iplot(fig, validate = False)
