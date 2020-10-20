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
black = [100/255, 100/255, 100/255]


#%%
path = '../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120_ADCT2w/'

eightwTumVolRec = np.loadtxt(path + "8wTumVolNormRec.res")
eightwTumVolNoRec = np.loadtxt(path + "8wTumVolNormNoRec.res")

stats.mannwhitneyu(eightwTumVolRec, eightwTumVolNoRec, alternative = 'less')



#%%

path = '../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120_ADCT2w/';
path20x3 = '../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120_ADCT2w_20x3/';

rec = [9, 10, 18, 39, 41, 47, 62, 69, 72];
nRec = 0;
nNoRec = 0;
tumVol = np.zeros([361, 2, 5])
tumVolRec = np.zeros([361, 2, 9])
tumVolNoRec = np.zeros([361, 2, 67])
for i in range(76) :
    if (i in rec) :     
        for k in range(5) :
            tumVol[:, :, k] = np.loadtxt(path20x3 + 'rep' + str(k) + '/tumVol_' +
                   str(i + 1) + '.res')
#        tumVol[:, 1, :] = tumVol[:, 1, :] / tumVol[0, 1, :]
        tumVolRec[:, :, nRec] = np.mean(tumVol, axis = 2)
        nRec = nRec + 1
    
    else :
        for k in range(5) :
            tumVol[:, :, k] = np.loadtxt(path + 'rep' + str(k) + '/tumVol_' +
                   str(i + 1) + '.res')
#        tumVol[:, 1, :] = tumVol[:, 1, :] / tumVol[0, 1, :]
        tumVolNoRec[:, :, nNoRec] = np.mean(tumVol, axis = 2)
        nNoRec = nNoRec + 1
        
meanTumVolRec = np.mean(tumVolRec, axis = 2);
stdTumVolRec = np.std(tumVolRec, axis = 2);
meanTumVolNoRec = np.mean(tumVolNoRec, axis = 2);
stdTumVolNoRec = np.std(tumVolNoRec, axis = 2);

meanTumVolRec[:, 0] = meanTumVolRec[:, 0] / (24 * 7)
meanTumVolNoRec[:, 0] = meanTumVolNoRec[:, 0] / (24 * 7)

#%%
fig, ax = plt.subplots() 

ax.plot(meanTumVolRec[:225, 0], meanTumVolRec[:225, 1], color = red, linewidth = 6)
ax.fill_between(meanTumVolRec[:225, 0],
                meanTumVolRec[:225, 1] - stdTumVolRec[:225, 1],
                meanTumVolRec[:225, 1] + stdTumVolRec[:225, 1], color = red, alpha = 0.1)
ax.plot(meanTumVolNoRec[:225, 0], meanTumVolNoRec[:225, 1], color = green, linewidth = 6)
ax.fill_between(meanTumVolNoRec[:225, 0],
                meanTumVolNoRec[:225, 1] - stdTumVolNoRec[:225, 1],
                meanTumVolNoRec[:225, 1] + stdTumVolNoRec[:225, 1], color = green, alpha = 0.1)

pvalue = np.ones(225)
for i in range(3, 225) :
    _, pvalue[i] = stats.mannwhitneyu(tumVolRec[i, 1, :].ravel(),
                                     tumVolNoRec[i, 1, :].ravel(),
                                     alternative = 'greater')

#ax.plot(meanTumVolRec[:, 0], pvalue, color = gray, linewidth = 6)
#ax.plot(meanTumVolRec[:, 0], pvalue < 0.001, color = gray, linewidth = 6)

ax.fill_between(meanTumVolRec[:225, 0], np.zeros(225), 1.1 * (pvalue <= 0.01),
                color = gray, linewidth = 0, alpha = 0.1)
ax.text(3.5, 1, 'p $\leq$ 0.001')

ax.plot([8, 8], [0, 1.1], '--', color = black, linewidth = 6)
ax.text(6.5, 1, 'Prediction 4', bbox = dict(facecolor = greenBlue + [0.5]))

ax.set(xlabel = 't (weeks)', ylabel = 'Normalized n° of tumor cells',
       xlim = [0, 8.1], ylim = [0, 1.1])
ax.set_xticks(np.linspace(0, 8, 9))
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'],
          loc = 'lower left')


#%%
fig, ax = plt.subplots() 

ax.plot(meanTumVolRec[:, 0], meanTumVolRec[:, 1], color = red, linewidth = 6)
ax.fill_between(meanTumVolRec[:, 0],
                meanTumVolRec[:, 1] - stdTumVolRec[:, 1],
                meanTumVolRec[:, 1] + stdTumVolRec[:, 1], color = red, alpha = 0.1)
ax.plot(meanTumVolNoRec[:, 0], meanTumVolNoRec[:, 1], color = green, linewidth = 6)
ax.fill_between(meanTumVolNoRec[:, 0],
                meanTumVolNoRec[:, 1] - stdTumVolNoRec[:, 1],
                meanTumVolNoRec[:, 1] + stdTumVolNoRec[:, 1], color = green, alpha = 0.1)

pvalue = np.ones(361)
for i in range(3, 361) :
    _, pvalue[i] = stats.mannwhitneyu(tumVolRec[i, 1, :].ravel(),
                                     tumVolNoRec[i, 1, :].ravel())

#ax.plot(meanTumVolRec[:, 0], pvalue, color = gray, linewidth = 6)
#ax.plot(meanTumVolRec[:, 0], pvalue < 0.001, color = gray, linewidth = 6)

ax.fill_between(meanTumVolRec[:, 0], np.zeros(361), 1.1 * (pvalue <= 0.001),
                color = gray, linewidth = 0, alpha = 0.1)
ax.text(3.5, 1, 'p $\leq$ 0.001')

ax.plot([8, 8], [0, 1.1], '--', color = black, linewidth = 6)
ax.text(6.5, 1, 'Prediction 4', bbox = dict(facecolor = greenBlue + [0.5]))

ax.set(xlabel = 't (weeks)', ylabel = 'Normalized n° of tumor cells',
       xlim = [0, 8.1], ylim = [0, 4.1])
ax.set_xticks(np.linspace(0, 12, 13))
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'],
          loc = 'upper right')
#%%
fig, ax = plt.subplots() 
ax.plot(meanTumVolRec[19:, 0], pvalue[19:], color = gray, linewidth = 6)



