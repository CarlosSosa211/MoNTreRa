#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np


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
black = [100/255, 100/255, 100/255]

#%%
path = '../../Carlos/Results/Recurrence/tissue26/'

nRep = 5
N = 361
tumVol = np.zeros((N, 2, nRep))
pO2 = np.zeros((N, 3, nRep))

for k in range(nRep) :
    tumVol[:, :, k] = np.loadtxt(path + 'rep' + str(k) + '/tumVol.res')
    pO2[:, :, k] = np.loadtxt(path + 'rep' + str(k) + '/pO2Stat.res')

meanTumVol = np.mean(tumVol[:, 1, :], axis = 1)
stdTumVol = np.std(tumVol[:, 1, :], axis = 1)
meanNormTumVol = np.mean(tumVol[:, 1, :] / tumVol[0, 1, :], axis = 1)
stdNormTumVol = np.std(tumVol[:, 1, :] / tumVol[0, 1, :], axis = 1)
meanMedPO2 = np.mean(pO2[:, 1, :], axis = 1)
stdMedPO2 = np.std(pO2[:, 1, :], axis = 1)

t = tumVol[:, 0, 0]
tW = t / (24 * 7)

#%%
N = 361
fig, ax = plt.subplots() 

cellSize = 0.02 * 0.02
ax.plot(tW[:N], meanTumVol[:N] / cellSize, color = blue, linewidth = 6)
ax.fill_between(tW[:N], (meanTumVol[:N] - stdTumVol[:N]) / cellSize,
                (meanTumVol[:N] + stdTumVol[:N]) / cellSize, color = blue,
                alpha = 0.25)
ax.set(xlabel = 't (weeks)', ylabel = 'N° of tumor cells',
       xlim = [0, 12.1])
ax.set_xticks(np.linspace(0, 12, 13))

#%%
N = 361
fig, ax = plt.subplots() 

ax.plot(tW[:N], meanNormTumVol[:N], color = blue, linewidth = 6)
ax.fill_between(tW[:N], meanNormTumVol[:N] - stdNormTumVol[:N],
                meanNormTumVol[:N] + stdNormTumVol[:N], color = blue,
                alpha = 0.25)
ax.set(xlabel = 't (weeks)', ylabel = 'Normalized n° of tumor cells',
       xlim = [0, 12.1])
ax.set_xticks(np.linspace(0, 12, 13))

#%%
N = 361
fig, ax = plt.subplots() 

cellSize = 0.02 * 0.02
ax.plot(tW[1:N], meanMedPO2[1:N], color = green, linewidth = 6)
ax.fill_between(tW[1:N], meanMedPO2[1:N] - stdMedPO2[1:N],
                meanMedPO2[1:N] + stdMedPO2[1:N],
                color = green, alpha = 0.25)
ax.set(xlabel = 't (weeks)', ylabel = 'Median $pO_2$ (mmHg)',
       xlim = [0, 12.1], ylim = [0, None])
ax.set_xticks(np.linspace(0, 12, 13))