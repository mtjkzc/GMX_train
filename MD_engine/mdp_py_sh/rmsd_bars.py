#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 08:28:48 2024

@author: mkozic
"""

###-----------------###
###  .xvg  grapher  ###
###-----------------###

def xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2):
    
    import numpy as np
    from numpy import loadtxt
    import pandas as pd
    import matplotlib.pyplot as plt
    from pandas import Series
    import re

    
    path = path
    filename = filename
    decimals = decimals
    color = color
    size = size
    savefig = savefig
    pos1 = pos1
    pos2 = pos2
    
    
    file = str(str(path) + str(filename))
    xvg_file= open(file)
    dataframe = pd.DataFrame(xvg_file)
    
    
    # DataFrame 1: Keep rows that start with "@"
    df1 = dataframe[dataframe.iloc[:,0].str.startswith('@')].copy()
    
    # DataFrame 2: Keep rows that start with "#"
    df2 = dataframe[dataframe.iloc[:,0].str.startswith('#')].copy()
    
    # DataFrame 3: Delete rows that contain "@" or "#"
    df3 = dataframe[~dataframe.iloc[:,0].str.contains('@|#')].copy()
    
    
    data = df3[0].str.split(expand=True)
    data=np.array(data, dtype="float")
    data=data.round(decimals)
    x = data[:,0]
    y = data[:,1]
    y = float(y[0])
    
    
    # Function to extract values for a given label
    def extractor(df, label):
        pattern = fr'{label}\s*"([^"]+)"'
        extracted_values = []
    
        for text in df[0]:
            match = re.search(pattern, text)
            if match:
                extracted_values.append(match.group(1))
            else:
                extracted_values.append(None)
                
        extracted_values = str(list(filter(None, extracted_values)))
        extracted_values = extracted_values.replace("[","")
        extracted_values = extracted_values.replace("]","")
        extracted_values = extracted_values.replace("'","") 
        return extracted_values
    
    
    title = extractor(df1, "title")
    label = extractor(df1, "legend")
    xlabel = extractor(df1, "xaxis  label")
    ylabel = extractor(df1, "yaxis  label")
    legend_loc = "upper left"
    
    out = [title, label, xlabel, ylabel, x, y, color, size, savefig, pos1, pos2]
    
    return out

###############################################################################
# START START START START START START START START START START START START STA #
###############################################################################

### Define the plot grid ###

import matplotlib.pyplot as plt
import numpy as np

path = "XnameX_MD_graphs/"
size = 800
linewidth = 0.4
savefig = "XnameX_fig6_rmsd_bars.png"
pos1 = pos2 = 0
color = "black"   
decimals = 8

A = "XnameX_rmsd_m0_m1_.xvg"
B = "XnameX_rmsd_m0_m2_.xvg"
C = "XnameX_rmsd_m1_m2_.xvg"

A_out = xvgreader(path, A, decimals, color, size, savefig, pos1, pos2)
B_out = xvgreader(path, B, decimals, color, size, savefig, pos1, pos2)
C_out = xvgreader(path, C, decimals, color, size, savefig, pos1, pos2)

A_value = A_out[5]
B_value = B_out[5]
C_value = C_out[5]

a = A_out[0] # title
b = A_out[1] # label
c = A_out[2] # xlabel
d = A_out[3] # ylabel

values = np.array([A_value, B_value, C_value])
labels = ["m0 vs m1", "m0 vs m2", "m1 vs m2"]

fig, ax = plt.subplots()

plt.bar(labels, values, color="black")

# function to add value labels
def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i, y[i], y[i], ha = 'center', color = "gray")
        
addlabels(x=labels, y=values)        

ax.set_ylabel(d)
ax.set_xlabel("Models compared")
ax.set_title("RMSD between AlphaFold structures")

plt.savefig(savefig, dpi=size)
#plt.show()

np.savetxt('XnameX_models_rmsd.csv', values, delimiter=',', header=','.join(labels), comments='')

###############################################################################
# END END END END END END END END END END END END END END END END END END END #
###############################################################################

