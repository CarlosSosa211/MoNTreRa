#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
import seaborn as sns
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
tcolor = [blue, blue, blue, blue, orange, blue, blue, orange, blue, redBlue]
rsc = [-1.31, -1.25, -1.25, -0.87, 0.80, -0.71, -0.48, 0.40, -0.26, -0.005]

fig, ax = plt.subplots()
ax = plt.bar(np.arange(len(rsc)), rsc, color = tcolor, 
             tick_label = ['$d_{vasc}$', r'$\alpha$', '$pO_2^{cap}$', '$D$', 
                           '$V_{max}$', r'$\beta$', '$K_{m}$', '$m$', '$t_{p}$',
                           '$t_{r}$'] ) 
plt.ylabel('RSC')
plt.grid(True, axis = 'y')
