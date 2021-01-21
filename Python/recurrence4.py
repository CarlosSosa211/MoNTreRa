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
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from scipy.fft import fft, ifft

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
path = '../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120/'
path28x25 = '../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120_28x2.5/'
path20x3 = '../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120_20x3/'

nTissues = 76
nRep = 10
N = 361
#rec = [9, 10, 18, 39, 41, 47, 62, 69, 72]
rec = [9, 10, 18, 25, 26, 31, 37, 39, 41, 47, 51, 57, 59, 62, 66, 69, 72]
nRec = len(rec)
nNoRec = nTissues - len(rec)
iRec = 0
iNoRec = 0

tumVol = np.zeros((N, 2, nRep))
tumVolRec = np.zeros((N, 2, nRep, nRec))
tumVolNoRec = np.zeros((N, 2, nRep, nNoRec))
for i in range(nTissues) :
    if (i in rec) :     
        for k in range(nRep) :
            tumVol[:, :, k] = np.loadtxt(path + 'rep' + str(k) + '/tumVol_' +
                   str(i + 1) + '.res')
        tumVol[:, 1, :] = tumVol[:, 1, :] / tumVol[0, 1, :]
        tumVolRec[:, :, :, iRec] = tumVol 
        iRec = iRec + 1
    
    else :
        for k in range(nRep) :
            tumVol[:, :, k] = np.loadtxt(path + 'rep' + str(k) + '/tumVol_' +
                   str(i + 1) + '.res')
        tumVol[:, 1, :] = tumVol[:, 1, :] / tumVol[0, 1, :]
        tumVolNoRec[:, :, :, iNoRec] = tumVol
        iNoRec = iNoRec + 1
        
meanTumVolRec = np.mean(tumVolRec, axis = 2)
stdTumVolRec = np.std(tumVolRec, axis = 2)
meanTumVolNoRec = np.mean(tumVolNoRec, axis = 2)
stdTumVolNoRec = np.std(tumVolNoRec, axis = 2)

meanTumVolRec[:, 0, :] = meanTumVolRec[:, 0, :] / (24 * 7)
meanTumVolNoRec[:, 0, :] = meanTumVolNoRec[:, 0, :] / (24 * 7)

w = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])
t = w * 28 #eight w

wTumVolRec = tumVolRec[t, 1, :, :]
wTumVolNoRec = tumVolNoRec[t, 1, :, :]

meanWTumVolRec = meanTumVolRec[t, 1, :]
stdWTumVolRec = stdTumVolRec[t, 1, :]

meanWTumVolNoRec = meanTumVolNoRec[t, 1, :]
stdWTumVolNoRec = stdTumVolNoRec[t, 1, :]

wIntTumVolRec = np.zeros((len(t), nRep, nRec))
wIntTumVolNoRec = np.zeros((len(t), nRep, nNoRec))

for i in range(len(t)) :
    wIntTumVolRec[i, :, :] = np.trapz(tumVolRec[:t[i], 1, :, :],
                 tumVolRec[:t[i], 0, :, :], axis = 0)
    wIntTumVolNoRec[i, :, :] = np.trapz(tumVolNoRec[:t[i], 1, :, :],
                   tumVolNoRec[:t[i], 0, :, :], axis = 0)

meanWIntTumVolRec = np.mean(wIntTumVolRec, axis = 1)
stdWIntTumVolRec = np.std(wIntTumVolRec, axis = 1)

meanWIntTumVolNoRec = np.mean(wIntTumVolNoRec, axis = 1)
stdWIntTumVolNoRec = np.std(wIntTumVolNoRec, axis = 1)


#%%
i = 2
ind = np.argsort(meanWTumVolRec[i, :])
meanWTumVolRec = meanWTumVolRec[:, ind]
stdWTumVolRec = stdWTumVolRec[:, ind]

ind = np.argsort(meanWTumVolNoRec[i, :])
meanWTumVolNoRec = meanWTumVolNoRec[:, ind]
stdWTumVolNoRec = stdWTumVolNoRec[:, ind]

ind = np.argsort(meanWIntTumVolRec[i, :])
meanWIntTumVolRec = meanWIntTumVolRec[:, ind]
stdWIntTumVolRec = stdWIntTumVolRec[:, ind]

ind = np.argsort(meanWIntTumVolNoRec[i, :])
meanWIntTumVolNoRec = meanWIntTumVolNoRec[:, ind]
stdWIntTumVolNoRec = stdWIntTumVolNoRec[:, ind]

#%%
fig, ax = plt.subplots() 
for k in range(5) :
    ax.scatter(np.linspace(0, nNoRec - 1, nNoRec), wTumVolNoRec[i, k, :],
               color = green)
    ax.scatter(np.linspace(nNoRec, 75, nRec), wTumVolRec[i, k, :],
               color = red)

#%%
for k in range(5, nRep) :
    ax.scatter(np.linspace(0, nNoRec - 1, nNoRec), wTumVolNoRec[i, k, :],
               color = greenBlue)
    ax.scatter(np.linspace(nNoRec, 75, nRec), wTumVolRec[i, k, :],
               color = redOrange)
    
#%%
fig, ax = plt.subplots() 
for k in range(nRep) :
    ax.scatter(np.linspace(0, nNoRec - 1, nNoRec), wIntTumVolNoRec[i, k, :],
               color = green)
    ax.scatter(np.linspace(nNoRec, 75, nRec), wIntTumVolRec[i, k, :],
               color = red)

#%%
fig, ax = plt.subplots() 
ax.errorbar(np.linspace(0, nNoRec - 1, nNoRec), meanWTumVolNoRec[i, :],
            stdWTumVolNoRec[i, :], marker = "o", fmt = ' ', color = green,
            ecolor = gray, elinewidth = 4)

ax.errorbar(np.linspace(nNoRec, 75, nRec), meanWTumVolRec[i, :],
            stdWTumVolRec[i, :], marker = "o", fmt = ' ', color = red,
            ecolor = gray, elinewidth = 4)
ax.grid()


#%%
fig, ax = plt.subplots() 
ax.errorbar(np.linspace(0, nNoRec - 1, nNoRec), meanWIntTumVolNoRec[i, :],
            stdWIntTumVolNoRec[i, :], marker = "o", fmt = ' ', color = green,
            ecolor = gray, elinewidth = 4)

ax.errorbar(np.linspace(nNoRec, 75, nRec), meanWIntTumVolRec[i, :],
            stdWIntTumVolRec[i, :], marker = "o", fmt = ' ', color = red,
            ecolor = gray, elinewidth = 4)
ax.grid()

#%%
fig, ax = plt.subplots() 
for k in range(5):
    ax.plot(tumVolNoRec[:, 0, k, 1], tumVolNoRec[:, 1, k, 1], linewidth = 6)


#%%
fig, ax = plt.subplots() 
ax.scatter(np.linspace(0, nNoRec - 1, nNoRec), meanWTumVolNoRec[0, :],
           color = greenBlue)
ax.scatter(np.linspace(nNoRec, 75, nRec), meanWTumVolRec[0, :],
           color = lightRed)

ax.scatter(np.linspace(0, nNoRec - 1, nNoRec), meanWTumVolNoRec[1, :],
           color = green)
ax.scatter(np.linspace(nNoRec, 75, nRec), meanWTumVolRec[1, :],
           color = red)

ax.scatter(np.linspace(0, nNoRec - 1, nNoRec), meanWTumVolNoRec[2, :],
           color = lightGreen)
ax.scatter(np.linspace(nNoRec, 75, nRec), meanWTumVolRec[2, :],
           color = redOrange)

ax.grid()