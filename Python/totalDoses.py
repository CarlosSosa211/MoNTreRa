#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
import seaborn as sns
from statannot import add_stat_annotation
from scipy.optimize import curve_fit
import scipy.stats as stats
    
#%%
plt.rcParams.update({'pdf.fonttype': 42})
plt.rcParams.update({'ps.fonttype': 42})
plt.rcParams.update({'font.size': 32})
green = [153/255, 255/255, 102/255]
darkGreen = [127/255, 207/255, 127/255]
lightGreen = [205/255, 255/255, 105/255]
blue = [127/255, 185/255, 225/255]
darkBlue = [150/255, 158/255, 196/255]
red = [255/255, 124/255, 128/255]
darkRed = [232/255, 136/255, 136/255]
lightRed = [255/255, 101/255, 152/255]
orange = [255/255, 174/255, 101/255]
gray = [173/255, 185/255, 202/255]
darkPurple = [160/255, 120/255, 196/255]
#lightPurple = [146/255, 141/255, 251/255]
lightPurple = [124/255, 169/255, 184/255]
redOrange = [255/255, 149/255, 114/255]
redPurple = [207/255, 122/255, 162/255]
redGreen = [230/255, 190/255, 117/255]
redBlue = [191/255, 155/255, 177/255]
orangePurple = [207/255, 147/255, 148/255]
greenBlue = [166/255, 220/255, 165/255]
brown = [168/255, 120/255, 110/255]

#%%
data = pd.read_csv('../../Carlos/Results/Hypofractionation/summary.csv')
data = data.rename(columns = {'fraction' : 'Fraction (Gy)'})
#%%
tcolor = [green, orange, red]

g = sns.catplot(x = "alphaBeta", y ="totalDose",
                hue = "fraction", col = "percentage",
                data = data, kind = "bar", palette = tcolor, col_wrap = 3);

#%%
fig, ax = plt.subplots()

tcolor = [green, orange, red]

perc = 99.9
dataPerc = data.loc[data['percentage'] == perc]
ax = plt.plot([-0.5, 5.5], [57, 57], color = gray, linewidth = 6,
              linestyle = 'dashed')
ax = plt.plot([-0.5, 5.5], [80, 80], color = gray, linewidth = 6,
              linestyle = 'dashed')
#ax = plt.scatter([3], [120], marker = '*', color = gray, s = 600)
#ax = plt.scatter([2], [200], marker = '*', color = gray, s = 600)
ax = plt.scatter([0, 1], [200, 200], marker = '*', color = gray, s = 600)

ax = sns.barplot(x = 'alphaBeta', y = 'totalDose', hue = 'Fraction (Gy)',
                 data = dataPerc, palette = tcolor, 
                 order = ['valdagni', 'walsh', 'wang', 'miralbell1', 'brenner',
                          'miralbell2']);
ax.legend([],[], frameon = False)
ax.set(xlabel = '', ylabel = 'Total dose (Gy)')
ax.set_xticklabels(['Valdagni\n$et$ $al$.', 'Walsh\n$et$ $al$.',
                    'Wang\n$et$ $al$.',
                    'Miralbell\n$et$ $al$. I', 'Brenner\n$et$ $al$.',
                    'Miralbell\n$et$ $al$. II'], ha = 'center')

