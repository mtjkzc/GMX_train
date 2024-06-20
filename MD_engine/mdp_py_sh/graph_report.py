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
    
    try:
        if "rmsf" in filename:
            i = 0
            while i < len(x):
                line = x[i]
                line2 = x[i+1]
                if line < line2:
                    i = i + 1
                elif line > line2:
                    x[i+1] = line + 1
                    i = i + 1
    except IndexError:
        pass
        
    
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

# MIN # 
import matplotlib.pyplot as plt
fig, ax = plt.subplots(3)

path = "XpdbnameX_MD_graphs/"
size = 800
linewidth = 0.4
savefig = "XpdbnameX_fig1_min.png"

###############################################################################
#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA#
###############################################################################

# MIN PDB
# Graph 0,0
filename = "XpdbnameX_m0_min_potential.xvg"
color = "black"   
decimals = 6
pos1 = 0
pos2 = 0

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j].plot(e,f, color = g, label = b, linewidth = linewidth)

ax[j].set_title(filename, fontsize = 8)
ax[j].set_xlabel(c)
ax[j].set_ylabel (d)
ax[j].xaxis.label.set_size(8)
ax[j].yaxis.label.set_size(8)
ax[j].tick_params(axis='both', which='major', labelsize=8)
ax[j].tick_params(axis='both', which='minor', labelsize=8)

##legend= plt.legend(fancybox=True, shadow=True)

fig.tight_layout()

#plt.savefig(savefig, dpi=size)

###############################################################################
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB#
###############################################################################

#MIN AF
# Graph 0,2
filename = "XpdbnameX_m1_min_potential.xvg"
color = "black"   
decimals = 6
pos1 = 1
pos2 = 1

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j].set_xlabel(c)
ax[j].set_ylabel (d)
ax[j].xaxis.label.set_size(8)
ax[j].yaxis.label.set_size(8)
ax[j].tick_params(axis='both', which='major', labelsize=8)
ax[j].tick_params(axis='both', which='minor', labelsize=8)

##legend= plt.legend(fancybox=True, shadow=True)

fig.tight_layout()

#plt.savefig(savefig, dpi=size)

###############################################################################
#CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC#
###############################################################################

#MIN m2
# Graph 0,3
filename = "XpdbnameX_m2_min_potential.xvg"
color = "black"   
decimals = 6
pos1 = 2
pos2 = 2

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j].set_xlabel(c)
ax[j].set_ylabel (d)
ax[j].xaxis.label.set_size(8)
ax[j].yaxis.label.set_size(8)
ax[j].tick_params(axis='both', which='major', labelsize=8)
ax[j].tick_params(axis='both', which='minor', labelsize=8)

##legend= plt.legend(fancybox=True, shadow=True)

fig.tight_layout()

plt.savefig(savefig, dpi=size)
plt.close

###############################################################################
#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD#
###############################################################################

# NVT # 
import matplotlib.pyplot as plt
fig, ax = plt.subplots(3,3)

path = "XpdbnameX_MD_graphs/"
size = 800
savefig = "XpdbnameX_fig2_nvt.png"

###############################################################################
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE#
###############################################################################

# NVT m0 T
# Graph 0,0
filename = "XpdbnameX_m0_nvt_temperature.xvg"
color = "black"   
decimals = 6
pos1 = 0
pos2 = 0

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)

ax[j,k].set_title(filename, fontsize = 8)
#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF#
###############################################################################
# NVT_m1_T
# Graph 0,1
filename = "XpdbnameX_m1_nvt_temperature.xvg"
color = "black"   
decimals = 6
pos1 = 0
pos2 = 1

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG#
###############################################################################

# NVT m2 T
# Graph 0,2
filename = "XpdbnameX_m2_nvt_temperature.xvg"
color = "black"   
decimals = 6
pos1 = 0
pos2 = 2

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
###############################################################################

# NVT m0 p
# Graph 1,0
filename = "XpdbnameX_m0_nvt_pressure.xvg"
color = "black"   
decimals = 6
pos1 = 1
pos2 = 0

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII#
###############################################################################
# NVT_m1_p
# Graph 1,1
filename = "XpdbnameX_m1_nvt_pressure.xvg"
color = "black"   
decimals = 6
pos1 = 1
pos2 = 1

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ#
###############################################################################
 
# NVT m0 p
# Graph 1,2
filename = "XpdbnameX_m2_nvt_pressure.xvg"
color = "black"   
decimals = 6
pos1 = 1
pos2 = 2

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK#
###############################################################################
# NVT m0 E
# Graph 2,0
filename = "XpdbnameX_m0_nvt_potential.xvg"
color = "black"   
decimals = 6
pos1 =2
pos2 = 0

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#R#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL#
###############################################################################
# NVT_m1_E
# Graph 2,1
filename = "XpdbnameX_m1_nvt_potential.xvg"
color = "black"   
decimals = 6
pos1 = 2
pos2 = 1

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM#
###############################################################################
# NVT m0 E
# Graph 2,2
filename = "XpdbnameX_m2_nvt_potential.xvg"
color = "black"   
decimals = 6
pos1 = 2
pos2 = 2

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN#
###############################################################################


# npt # 
import matplotlib.pyplot as plt
fig, ax = plt.subplots(3,3)

path = "XpdbnameX_MD_graphs/"
size = 800
savefig = "XpdbnameX_fig3_npt.png"

###############################################################################
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE#
###############################################################################

# npt m0 T
# Graph 0,0
filename = "XpdbnameX_m0_npt_temperature.xvg"
color = "black"   
decimals = 6
pos1 = 0
pos2 = 0

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)

ax[j,k].set_title(filename, fontsize = 8)
#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF#
###############################################################################
# npt_m1_T
# Graph 0,1
filename = "XpdbnameX_m1_npt_temperature.xvg"
color = "black"   
decimals = 6
pos1 = 0
pos2 = 1

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG#
###############################################################################

# npt m2 T
# Graph 0,2
filename = "XpdbnameX_m2_npt_temperature.xvg"
color = "black"   
decimals = 6
pos1 = 0
pos2 = 2

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
###############################################################################

# npt m0 p
# Graph 1,0
filename = "XpdbnameX_m0_npt_pressure.xvg"
color = "black"   
decimals = 6
pos1 = 1
pos2 = 0

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII#
###############################################################################
# npt_m1_p
# Graph 1,1
filename = "XpdbnameX_m1_npt_pressure.xvg"
color = "black"   
decimals = 6
pos1 = 1
pos2 = 1

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ#
###############################################################################
 
# npt m0 p
# Graph 1,2
filename = "XpdbnameX_m2_npt_pressure.xvg"
color = "black"   
decimals = 6
pos1 = 1
pos2 = 2

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK#
###############################################################################
# npt m0 E
# Graph 2,0
filename = "XpdbnameX_m0_npt_potential.xvg"
color = "black"   
decimals = 6
pos1 =2
pos2 = 0

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#R#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL#
###############################################################################
# npt_m1_E
# Graph 2,1
filename = "XpdbnameX_m1_npt_potential.xvg"
color = "black"   
decimals = 6
pos1 = 2
pos2 = 1

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM#
###############################################################################
# npt m0 E
# Graph 2,2
filename = "XpdbnameX_m2_npt_potential.xvg"
color = "black"   
decimals = 6
pos1 = 2
pos2 = 2

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO#
###############################################################################


# md # 
import matplotlib.pyplot as plt
fig, ax = plt.subplots(3,3)

path = "XpdbnameX_MD_graphs/"
size = 800
savefig = "XpdbnameX_fig4_md.png"

###############################################################################
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE#
###############################################################################

# md m0 T
# Graph 0,0
filename = "XpdbnameX_m0_md_temperature.xvg"
color = "black"   
decimals = 6
pos1 = 0
pos2 = 0

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)

ax[j,k].set_title(filename, fontsize = 8)
#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF#
###############################################################################
# md_m1_T
# Graph 0,1
filename = "XpdbnameX_m1_md_temperature.xvg"
color = "black"   
decimals = 6
pos1 = 0
pos2 = 1

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG#
###############################################################################

# md m2 T
# Graph 0,2
filename = "XpdbnameX_m2_md_temperature.xvg"
color = "black"   
decimals = 6
pos1 = 0
pos2 = 2

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
###############################################################################

# md m0 p
# Graph 1,0
filename = "XpdbnameX_m0_md_pressure.xvg"
color = "black"   
decimals = 6
pos1 = 1
pos2 = 0

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII#
###############################################################################
# md_m1_p
# Graph 1,1
filename = "XpdbnameX_m1_md_pressure.xvg"
color = "black"   
decimals = 6
pos1 = 1
pos2 = 1

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ#
###############################################################################
 
# md m0 p
# Graph 1,2
filename = "XpdbnameX_m2_md_pressure.xvg"
color = "black"   
decimals = 6
pos1 = 1
pos2 = 2

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK#
###############################################################################
# md m0 E
# Graph 2,0
filename = "XpdbnameX_m0_md_potential.xvg"
color = "black"   
decimals = 6
pos1 =2
pos2 = 0

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#R#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL#
###############################################################################
# md_m1_E
# Graph 2,1
filename = "XpdbnameX_m1_md_potential.xvg"
color = "black"   
decimals = 6
pos1 = 2
pos2 = 1

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM#
###############################################################################
# md m0 E
# Graph 2,2
filename = "XpdbnameX_m2_md_potential.xvg"
color = "black"   
decimals = 6
pos1 = 2
pos2 = 2

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP#
###############################################################################


# md results # 
import matplotlib.pyplot as plt
fig, ax = plt.subplots(3,3)

path = "XpdbnameX_MD_graphs/"
size = 800
savefig = "XpdbnameX_fig5_md_results.png"

###############################################################################
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE#
###############################################################################

# md m0 rmsd
# Graph 0,0
filename = "XpdbnameX_m0_rmsd.xvg"
color = "black"   
decimals = 6
pos1 = 0
pos2 = 0

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)

ax[j,k].set_title(filename, fontsize = 8)
#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF#
###############################################################################
# md_m1_rmsd
# Graph 0,1
filename = "XpdbnameX_m1_rmsd.xvg"
color = "black"   
decimals = 6
pos1 = 0
pos2 = 1

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG#
###############################################################################

# md m2 rmsd
# Graph 0,2
filename = "XpdbnameX_m2_rmsd.xvg"
color = "black"   
decimals = 6
pos1 = 0
pos2 = 2

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
###############################################################################

# md m0 rmsf
# Graph 1,0
filename = "XpdbnameX_m0_rmsf.xvg"
color = "black"   
decimals = 6
pos1 = 1
pos2 = 0

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII#
###############################################################################
# md_m1_rmsf
# Graph 1,1
filename = "XpdbnameX_m1_rmsf.xvg"
color = "black"   
decimals = 6
pos1 = 1
pos2 = 1

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ#
###############################################################################
 
# md m0 rmsf
# Graph 1,2
filename = "XpdbnameX_m2_rmsf.xvg"
color = "black"   
decimals = 6
pos1 = 1
pos2 = 2

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK#
###############################################################################
# md m0 gyr
# Graph 2,0
filename = "XpdbnameX_m0_gyr.xvg"
color = "black"   
decimals = 6
pos1 =2
pos2 = 0

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#R#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL#
###############################################################################
# md_m1_gyr
# Graph 2,1
filename = "XpdbnameX_m1_gyr.xvg"
color = "black"   
decimals = 6
pos1 = 2
pos2 = 1

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

#plt.savefig(savefig, dpi=size)
fig.tight_layout()
###############################################################################
#MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM#
###############################################################################
# md m0 gyr
# Graph 2,2
filename = "XpdbnameX_m2_gyr.xvg"
color = "black"   
decimals = 6
pos1 = 2
pos2 = 2

out = xvgreader(path, filename, decimals, color, size, savefig, pos1, pos2)
a = out[0] # title
b = out[1] # label
c = out[2] # xlabel
d = out[3] # ylabel
e = out[4] # x
f = out[5] # y
g = out[6] # color
h = out[7] # size
i = out[8] # savefig
j = out[9] # pos1
k = out[10] # pos2
ax[j,k].plot(e,f, color = g, label = b, linewidth = linewidth)
ax[j,k].set_title(filename, fontsize = 8)

#plt.title(a)
ax[j,k].set_xlabel(c)
ax[j,k].set_ylabel (d)
ax[j,k].xaxis.label.set_size(6)
ax[j,k].yaxis.label.set_size(6)
ax[j,k].tick_params(axis='both', which='major', labelsize=5)
ax[j,k].tick_params(axis='both', which='minor', labelsize=5)

#legend= plt.legend(fancybox=True, shadow=True)

plt.savefig(savefig, dpi=size)
fig.tight_layout()

###############################################################################
#QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ#
###############################################################################


###############################################################################
# END END END END END END END END END END END END END END END END END END END #
###############################################################################

