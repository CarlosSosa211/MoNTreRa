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
plt.close('all')
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
path = '../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120/'

eightwTumVolRec = np.loadtxt(path + "8wTumVolNormRec.res")
eightwTumVolNoRec = np.loadtxt(path + "8wTumVolNormNoRec.res")

stats.mannwhitneyu(eightwTumVolRec, eightwTumVolNoRec, alternative = 'less')

#%%
path = '../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120/'
path28x25 = '../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120_28x2.5/'
path20x3 = '../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120_20x3/'

nTissues = 76
nRep = 5
N = 361
dx = 6
#rec = [9, 10, 18, 39, 41, 47, 62, 69, 72]
tRec = [9, 10, 18, 25, 26, 31, 37, 39, 41, 47, 51, 57, 59, 62, 66, 69, 72]
nRec = len(tRec)
nNoRec = nTissues - len(tRec)
iRec = 0
iNoRec = 0
rec = np.zeros(nTissues)
rec[tRec] = 1 
tNoRec = np.where(rec == 0)[0].tolist()
tumVolK = np.zeros((N, 2, nRep))
tumVol = np.zeros((N, nTissues))

k = 3
kern = np.ones(2 * k + 1) / (2 * k + 1)

for i in range(nTissues) :   
    for k in range(nRep) :
        tumVolK[:, :, k] = np.loadtxt(path + 'rep' + str(k) + '/tumVol_' +
              str(i + 1) + '.res')
#        tumVolK[:, 1, :] = tumVolK[:, 1, :] / tumVolK[0, 1, :]
    tumVol[:, i] = np.mean(tumVolK[:, 1, :], axis = 1)
#   tumVol[:, 1, i] = np.convolve(tumVol[:, 1, i], kern, mode = 'same')

t = tumVolK[:, 0, 0]
tW = t / (24 * 7)

meanTumVolRec = np.mean(tumVol[:, tRec], axis = 1)
stdTumVolRec = np.std(tumVol[:, tRec], axis = 1)

meanTumVolNoRec = np.mean(tumVol[:, tNoRec], axis = 1)
stdTumVolNoRec = np.std(tumVol[:, tNoRec], axis = 1)

wn = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
tn = wn * 28 #eight w

wTumVol = tumVol[tn, :]

meanWTumVolRec = meanTumVolRec[tn]
stdWTumVolRec = stdTumVolRec[tn]

meanWTumVolNoRec = meanTumVolNoRec[tn]
stdWTumVolNoRec = stdTumVolNoRec[tn]

wIntTumVol = np.zeros((len(tn), nTissues))

for i in range(len(tn)) :
    wIntTumVol[i, :] = np.trapz(tumVol[:tn[i], :], dx = dx, axis = 0)

meanWIntTumVolRec = np.mean(wIntTumVol[:, tRec], axis = 1)
stdWIntTumVolRec = np.std(wIntTumVol[:, tRec], axis = 1)

meanWIntTumVolNoRec = np.mean(wIntTumVol[:, tNoRec], axis = 1)
stdWIntTumVolNoRec = np.std(wIntTumVol[:, tNoRec], axis = 1)

pvalueTumVol = np.ones(N)
for i in range(3, N) :
    _, pvalueTumVol[i] = stats.mannwhitneyu(tumVol[i, tRec].ravel(),
                 tumVol[i, tNoRec].ravel())

#%%
fig, ax = plt.subplots() 

tcolor = [greenBlue, lightRed, green, red, lightGreen, redOrange]

wRec = np.tile(rec, len(wn)).reshape(nTissues * len(wn), 1)

twn = np.zeros((len(wn) * nTissues, 1))
for i in range(len(wn)) :
    twn[i * nTissues: (i + 1) * nTissues] = wn[i]

wTumVol_ = np.reshape(wTumVol, len(twn), order = 'C')

data = np.concatenate((wRec, twn, wTumVol_.reshape((len(wTumVol_), 1))),
                      axis = 1)
wTumVolDF = pd.DataFrame(data = data, columns = ['bio_rec', 'w',
                                                 'norm_tum_vol'])

ax = sns.boxplot(data = wTumVolDF, x = 'w', y = 'norm_tum_vol', hue = 'bio_rec',
                 orient = 'v', palette = tcolor)

#%%
i = 12
ind = np.argsort(wTumVol[i, :])

fig, ax = plt.subplots()
for j, ii in enumerate(ind) :
    if ii in tNoRec :
        ax.scatter(j, wTumVol[i, ii], color = green)
    else : 
        ax.scatter(j, wTumVol[i, ii], color = red)
ax.grid()

#%%
fig, ax = plt.subplots() 

tcolor = [greenBlue, lightRed, green, red, lightGreen, redOrange]

wRec = np.tile(rec, len(wn)).reshape(nTissues * len(wn), 1)

twn = np.zeros((len(wn) * nTissues, 1))
for i in range(len(wn)) :
    twn[i * nTissues: (i + 1) * nTissues] = wn[i]
    
wIntTumVol_ = np.reshape(wIntTumVol, len(twn), order = 'C')
data = np.concatenate((wRec, twn, wIntTumVol_.reshape((len(wIntTumVol_)), 1)),
                      axis = 1)
wIntTumVolDF = pd.DataFrame(data = data,
                            columns = ['bio_rec', 'w', 'int_norm_tum_vol'])

ax = sns.boxplot(data = wIntTumVolDF, x = 'w', y = 'int_norm_tum_vol',
                 hue = 'bio_rec', orient = 'v', palette = tcolor)

#%%
i = 5
ind = np.argsort(wIntTumVol[i, :])

fig, ax = plt.subplots()
for j, ii in enumerate(ind) :
    if ii in tNoRec :
        ax.scatter(j, wIntTumVol[i, ii], color = green)
    else : 
        ax.scatter(j, wIntTumVol[i, ii], color = red)
ax.grid()

#%%
N = 361
fig, ax = plt.subplots() 

pltNoRec = ax.plot(tW[:N], tumVol[:N, tNoRec], color = green, linewidth = 6,
                   alpha = 0.5, label = 'No biochemical recurrence')
plt.setp(pltNoRec[1:], label="_")
pltRec = ax.plot(tW[:N], tumVol[:N, tRec], color = red, linewidth = 6,
                 alpha = 0.5, label = 'Biochemical recurrence')
plt.setp(pltRec[1:], label="_")
ax.set(xlabel = 't (weeks)', ylabel = 'N° of tumor cells',
       xlim = [0, 12.1])
ax.set_xticks(np.linspace(0, 12, 13))
ax.legend(loc = 'upper right')
    
#%%
fig, ax = plt.subplots() 
N = 361

ax.plot(tW[:N], meanTumVolRec[:N], color = red, linewidth = 6)
ax.fill_between(tW[:N], meanTumVolRec[:N] - stdTumVolRec[:N],
                meanTumVolRec[:N] + stdTumVolRec[:N], color = red,
                alpha = 0.1)
ax.plot(tW[:N], meanTumVolNoRec[:N], color = green, linewidth = 6)
ax.fill_between(tW[:N], meanTumVolNoRec[:N] - stdTumVolNoRec[:N],
                meanTumVolNoRec[:N] + stdTumVolNoRec[:N], color = green,
                alpha = 0.1)

#ax.plot(meanTumVolRec[:, 0], pvalue, color = gray, linewidth = 6)


upLim = 1.1 * max(meanTumVolRec[:N] + stdTumVolRec[:N]) * np.ones(N)
lowLim = 0.9 * min(meanTumVolNoRec[:N] - stdTumVolNoRec[:N]) * np.ones(N)
ax.fill_between(tW[:N], lowLim, upLim, where = pvalueTumVol[:N] <= 0.001,
                color = gray, linewidth = 0, alpha = 0.5)
#ax.text(3.5, 1, 'p $\leq$ 0.001')
#
#ax.plot([8, 8], [0, 1.1], '--', color = black, linewidth = 6)
#ax.text(6.5, 1, 'Prediction 4', bbox = dict(facecolor = greenBlue + [0.5]))

ax.set(xlabel = 't (weeks)', ylabel = 'N° of tumor cells',
       xlim = [0, 12.1])
ax.set_xticks(np.linspace(0, 12, 13))
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'],
          loc = 'upper right')

#%%
fig, ax = plt.subplots() 
ax.plot(t[19:], pvalueTumVol[19:], color = gray, linewidth = 6)

#%%%
Nd = N - 1

dNoRec = np.zeros((Nd, nNoRec))
dRec = np.zeros((Nd, nRec))
dx = 6

k = 5
kern = np.ones(2 * k + 1) / (2 * k + 1)

diff = np.diff(tumVol, axis = 0) / dx
for i in range(nTissues) :
    diff[:, i] = np.convolve(diff[:, i], kern, mode = 'same')
        
wn = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
tn = wn * 28 #eight w

wDiff = diff[tn, :]

meanDiffNoRec = np.mean(diff[:, tNoRec], axis = 1)
stdDiffNoRec = np.std(diff[:, tNoRec], axis = 1)

meanDiffRec = np.mean(diff[:, tRec], axis = 1)
stdDiffRec = np.std(diff[:, tRec], axis = 1)

wIntDiff = np.zeros((len(tn), nTissues))

for i in range(len(tn)) :
    wIntDiff[i, :] = np.trapz(diff[:tn[i], :], dx = dx, axis = 0)

meanWIntDiffRec = np.mean(wIntDiff[:, tRec], axis = 1)
stdWIntDiffRec = np.std(wIntDiff[:, tRec], axis = 1)

meanWIntTumVolNoRec = np.mean(wIntDiff[:, tNoRec], axis = 1)
stdWIntTumVolNoRec = np.std(wIntDiff[:, tNoRec], axis = 1)

pvalueDiff = np.ones(Nd)
for i in range(3, Nd) :
    _, pvalueDiff[i] = stats.mannwhitneyu(diff[i, tRec].ravel(),
                 diff[i, tNoRec].ravel())

#%%
fig, ax = plt.subplots() 

tcolor = [greenBlue, lightRed, green, red, lightGreen, redOrange]

wRec = np.tile(rec, len(wn)).reshape(nTissues * len(wn), 1)

twn = np.zeros((len(wn) * nTissues, 1))
for i in range(len(wn)) :
    twn[i * nTissues: (i + 1) * nTissues] = wn[i]

wDiff_ = np.reshape(wDiff, len(twn), order = 'C')

data = np.concatenate((wRec, twn, wDiff_.reshape((len(wDiff_), 1))),
                      axis = 1)
wDiffDF = pd.DataFrame(data = data, columns = ['bio_rec', 'w',
                                               'diff_norm_tum_vol'])

ax = sns.boxplot(data = wDiffDF, x = 'w', y = 'diff_norm_tum_vol',
                 hue = 'bio_rec', orient = 'v', palette = tcolor)

#%%
i = -1
ind = np.argsort(wDiff[i, :])

fig, ax = plt.subplots()
for j, ii in enumerate(ind) :
    if ii in tNoRec :
        ax.scatter(j, wDiff[i, ii], color = green)
    else : 
        ax.scatter(j, wDiff[i, ii], color = red)
ax.grid()

#%%
fig, ax = plt.subplots() 

tcolor = [greenBlue, lightRed, green, red, lightGreen, redOrange]

wRec = np.tile(rec, len(wn)).reshape(nTissues * len(wn), 1)

twn = np.zeros((len(wn) * nTissues, 1))
for i in range(len(wn)) :
    twn[i * nTissues: (i + 1) * nTissues] = wn[i]

wIntDiff_ = np.reshape(wIntDiff, len(twn), order = 'C')

data = np.concatenate((wRec, twn, wIntDiff_.reshape((len(wIntDiff_), 1))),
                      axis = 1)
wIntDiffDF = pd.DataFrame(data = data, columns = ['bio_rec', 'w',
                                                  'int_diff_norm_tum_vol'])
ax = sns.boxplot(data = wIntDiffDF, x = 'w', y = 'int_diff_norm_tum_vol',
                 hue = 'bio_rec', orient = 'v', palette = tcolor)

#%%
i = -1
ind = np.argsort(wIntDiff[i, :])

fig, ax = plt.subplots() 
for j, ii in enumerate(ind) :
    if ii in tNoRec :
        ax.scatter(j, wIntDiff[i, ii], color = green)
    else : 
        ax.scatter(j, wIntDiff[i, ii], color = red)
ax.grid()

#%%
Nd = 360
fig, ax = plt.subplots() 

pltNoRec = ax.plot(tW[:Nd], diff[:Nd, tNoRec], color = green, linewidth = 6,
                   alpha = 0.5, label = 'No biochemical recurrence')
plt.setp(pltNoRec[1:], label="_")
pltRec = ax.plot(tW[:Nd], diff[:Nd, tRec], color = red, linewidth = 6,
                 alpha = 0.5, label = 'Biochemical recurrence')
plt.setp(pltRec[1:], label="_")
ax.set(xlabel = 't (weeks)', ylabel = 'Diff. of n° of tumor cells',
       xlim = [0, 12.1])
ax.set_xticks(np.linspace(0, 12, 13))
ax.legend(loc = 'upper left')

#%%
fig, ax = plt.subplots() 
Nd = 360

ax.plot(tW[:Nd], meanDiffRec[:Nd], color = red, linewidth = 6)
ax.fill_between(tW[:Nd], meanDiffRec[:Nd] - stdDiffRec[:Nd],
                meanDiffRec[:Nd] + stdDiffRec[:Nd], color = red, alpha = 0.1)
ax.plot(tW[:Nd], meanDiffNoRec[:Nd], color = green, linewidth = 6)
ax.fill_between(tW[:Nd], meanDiffNoRec[:Nd] - stdDiffNoRec[:Nd],
                meanDiffNoRec[:Nd] + stdDiffNoRec[:Nd], color = green,
                alpha = 0.1)

#ax.plot(meanTumVolRec[:, 0], pvalue, color = gray, linewidth = 6)

upLim = 1.1 * max(meanDiffRec[:Nd] + stdDiffRec[:Nd]) * np.ones(Nd)
lowLim = 0.9 * min(meanDiffNoRec[:Nd] - stdDiffNoRec[:Nd]) * np.ones(Nd)
ax.fill_between(tW[:Nd], lowLim, upLim, where = pvalueDiff[:Nd] <= 0.001,
                color = gray, linewidth = 0, alpha = 0.5)

ax.set(xlabel = 't (weeks)', ylabel = 'Diff. of n° of tumor cells',
       xlim = [0, 12.1])
ax.set_xticks(np.linspace(0, 12, 13))
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'],
          loc = 'upper left')
    
#%%
i = 3
ind = np.argsort(wIntTumVol[i, :])

fig, ax = plt.subplots()
for j, ii in enumerate(ind) :
    if ii in tNoRec :
        ax.scatter(wTumVol[i, ii], wIntTumVol[i, ii], color = green)
    else : 
        ax.scatter(wTumVol[i, ii], wIntTumVol[i, ii], color = red)
ax.grid()

#%%%
T = 6
Nd = N - 1

fTumVol = np.linspace(0.0, 1.0 / (2.0 * T), N // 2)
fDiff = np.linspace(0.0, 1.0 / (2.0 * T), Nd // 2)

fftTumVol = np.zeros((N, nTissues))
fftDiff = np.zeros((Nd, nTissues))

for i in range(nTissues):
    fftTumVol[:, i] = fft(tumVol[:N, i])
    fftDiff[:, i] = fft(diff[:Nd, i])
#    fftDiff[:, i] = fft(diff[:Nd,, i] + np.sin(2.0 * np.pi * 0.05 * t))
    
fftTumVol = 2.0 / N * np.abs(fftTumVol[0:N // 2, :])
fftDiff = 2.0 / Nd * np.abs(fftDiff[0:Nd // 2, :])
    
meanFftTumVolNoRec = np.mean(fftTumVol[:, tNoRec], axis = 1)
stdFftTumVolNoRec = np.std(fftTumVol[:, tNoRec], axis = 1)

meanFftTumVolRec = np.mean(fftTumVol[:, tRec], axis = 1)
stdFftTumVolRec = np.std(fftTumVol[:, tRec], axis = 1)

meanFftDiffNoRec = np.mean(fftDiff[:, tNoRec], axis = 1)
stdFftDiffNoRec = np.std(fftDiff[:, tNoRec], axis = 1)

meanFftDiffRec = np.mean(fftDiff[:, tRec], axis = 1)
stdFftDiffRec = np.std(fftDiff[:, tRec], axis = 1)

pvalueFftTumVol = np.ones(N // 2)
for i in range(N // 2) :
    _, pvalueFftTumVol[i] = stats.mannwhitneyu(fftTumVol[i, tRec].ravel(),
                      fftTumVol[i, tNoRec].ravel())

pvalueFftDiff = np.ones(Nd // 2)
for i in range(Nd // 2) :
    _, pvalueFftDiff[i] = stats.mannwhitneyu(fftDiff[i, tRec].ravel(),
                    fftDiff[i, tNoRec].ravel())

#%%
fig, ax = plt.subplots() 
Nf = 180
ax.plot(fTumVol[:Nf], meanFftTumVolRec[:Nf], color = red, linewidth = 6)
ax.fill_between(fTumVol[:Nf], meanFftTumVolRec[:Nf] - stdFftTumVolRec[:Nf],
                meanFftTumVolRec[:Nf] + stdFftTumVolRec[:Nf],
                color = red, alpha = 0.1)

ax.plot(fTumVol[:Nf], meanFftTumVolNoRec[:Nf], color = green, linewidth = 6)
ax.fill_between(fTumVol[:Nf], meanFftTumVolNoRec[:Nf] - stdFftTumVolNoRec[:Nf],
                meanFftTumVolNoRec[:Nf] + stdFftTumVolNoRec[:Nf],
                color = green, alpha = 0.1)

#ax.plot(meanTumVolRec[:, 0], pvalue, color = gray, linewidth = 6)

upLim = 1.1 * max(meanFftTumVolRec[:Nf] + stdFftTumVolRec[:Nf]) * np.ones(Nf)
lowLim = 0.9 * min(meanFftTumVolNoRec[:Nf] -
                   stdFftTumVolNoRec[:Nf]) * np.ones(Nf)
ax.fill_between(fDiff[:Nf], lowLim, upLim, where = pvalueFftTumVol[:Nf] <= 0.001,
                color = gray, linewidth = 0, alpha = 0.5)
    
ax.set(xlabel = 'f', ylabel = 'Fft. of normalized n° of tumor cells')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'],
          loc = 'upper right')

#%%
fig, ax = plt.subplots() 
Nf = 18
ax.plot(fDiff[:Nf], meanFftDiffRec[:Nf], color = red, linewidth = 6)
ax.fill_between(fDiff[:Nf], meanFftDiffRec[:Nf] - stdFftDiffRec[:Nf],
                meanFftDiffRec[:Nf] + stdFftDiffRec[:Nf],
                color = red, alpha = 0.1)

ax.plot(fDiff[:Nf], meanFftDiffNoRec[:Nf], color = green, linewidth = 6)
ax.fill_between(fDiff[:Nf], meanFftDiffNoRec[:Nf] - stdFftDiffNoRec[:Nf],
                meanFftDiffNoRec[:Nf] + stdFftDiffNoRec[:Nf],
                color = green, alpha = 0.1)

#ax.plot(meanDiffRec[:, 0], pvalue, color = gray, linewidth = 6)

upLim = 1.1 * max(meanFftDiffRec[:Nf] + stdFftDiffRec[:Nf]) * np.ones(Nf)
lowLim = 0.9 * min(meanFftDiffNoRec[:Nf] - stdFftDiffNoRec[:Nf]) * np.ones(Nf)
ax.fill_between(fDiff[:Nf], lowLim, upLim, where = pvalueDiff[:Nf] <= 0.001,
                color = gray, linewidth = 0, alpha = 0.5)

ax.set(xlabel = 'f', ylabel = 'Fft. of diff. of n° of tumor cells')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'],
          loc = 'upper right')
    
#%%
fig, ax = plt.subplots() 
Nf = 10
for i in tNoRec :
    ax.plot(fTumVol[:Nf], fftTumVol[:Nf, i], color = green, linewidth = 6,
            alpha = 0.5)
    
for i in tRec :
    ax.plot(fTumVol[:Nf], fftTumVol[:Nf, i], color = red, linewidth = 6,
            alpha = 0.5)    
    
ax.set(xlabel = 'f', ylabel = 'Fft. of normalized n° of tumor cells')
#ax.legend(['Biochemical recurrence', 'No biochemical recurrence'],
#          loc = 'upper right')

#%%
fig, ax = plt.subplots() 
Nf = 25
for i in tNoRec :
    ax.plot(fDiff[:Nf], fftDiff[:Nf, i], color = green, linewidth = 6,
            alpha = 0.5)
    
for i in tRec :
    ax.plot(fDiff[:Nf], fftDiff[:Nf, i], color = red, linewidth = 6,
            alpha = 0.5)    
    
ax.set(xlabel = 'f', ylabel = 'Fft. of diff. of normalized n° of tumor cells')
#ax.legend(['Biochemical recurrence', 'No biochemical recurrence'],
#          loc = 'upper right')
#%%
data = pd.read_csv('../../Carlos/Results/Recurrence/simp/rec_summary_frac.csv')
N = 1
#rec = [9, 10, 18, 39, 41, 47, 62, 69, 72]
rec = [9, 10, 18, 25, 26, 31, 37, 39, 41, 47, 51, 57, 59, 62, 66, 69, 72]
noRec = []
nRec = len(rec)
nNoRec = nTissues - len(rec)
for i in range(nTissues) :
    if i not in rec :
        noRec.append(i)
        
#tx = [data[['8w_int_tum_area']],
#      data[['8w_int_tum_area_20x3']]]

tx = [data[['8w_int_tum_area']],
      data[['8w_int_tum_area_28x2.5']],
      data[['8w_int_tum_area_20x3']]]

y = data[['bio_rec_6']].to_numpy().ravel() 

probas = np.zeros((nTissues, 2, len(tx)))
ypred = np.zeros((nTissues, len(tx)))
 
clf = LogisticRegression()
#clf = RandomForestClassifier()
#clf = MLPClassifier(activation = 'logistic', max_iter = 1000)

#threshold = 0.09
threshold = 0.25
 
for j in range(N) :
    x = tx[0].to_numpy()
    clf.fit(x, y)
    for i, x in enumerate(tx) :
        x = x.to_numpy()
        probas[:, :, i] += clf.predict_proba(x)
        ypred[:, i] += probas[:, 1, i] > threshold

probas /= N
ypred /= N 

ind = np.argsort(probas[:, 1, 0])
probas = probas[ind, :, :]
y = y[ind]
ypred = ypred[ind, :]

#medProbNoRec = np.median(probas[y == 0, 1, :], axis = 0)
#q1ProbNoRec = np.quantile(probas[y == 0, 1, :], 0.25, axis = 0)
#q3ProbNoRec = np.quantile(probas[y == 0, 1, :], 0.75, axis = 0)
#medProbRec = np.median(probas[y == 1, 1, :], axis = 0)

#%%
fig, ax = plt.subplots()
ax.fill_between([-1, nTissues], [0, 0], [threshold, threshold], color = green,
                alpha = 0.3)
ax.fill_between([-1, nTissues], [threshold, threshold], [1, 1], color = red,
                alpha = 0.3)
ax.scatter(np.linspace(0, nNoRec - 1, nNoRec), probas[y == 0, 1, 0],
           color = green, label = 'No recurrence (Standard fractionation)')
ax.scatter(np.linspace(nNoRec, 75, nRec), probas[y == 1, 1, 0], color = red,
           label = 'Recurrence (Standard fractionation)')
#ax.scatter(np.linspace(67, 75, 9), probas[rec, 1, 1], marker= 'v', color = green,
#           label = '28 x 2.5 Gy')
ax.scatter(np.linspace(nNoRec, 75, nRec), probas[y == 1, 1, 1], marker = 'v',
           color = red, label = 'Recurrence (28 x 2.5 Gy)')
#ax.scatter(np.linspace(67, 75, 9), probas[rec, 1, 2], marker = 's', color = green,
#           label = '20 x 3 Gy')
ax.scatter(np.linspace(nNoRec, 75, nRec), probas[y == 1, 1, 2], marker = 's',
           color = red, label = 'Recurrence (20 x 3 Gy)')

#ax.scatter(np.linspace(0, 66, 67), probas[noRec, 1, 1], color = blue)
#ax.fill_between([0, 67], [0, 0], [1, 1], color = green, alpha = 0.3)
#ax.fill_between([67, 75], [0, 0], [1, 1], color = red, alpha = 0.3)
ax.set(ylabel = 'Probability', xticks = [])
ax.legend(loc = 'upper left')

#%%
fig, ax = plt.subplots()
ax.fill_between([-1, nRec], [0, 0], [threshold, threshold], color = green,
                alpha = 0.3)
ax.fill_between([-1, nRec], [threshold, threshold], [1, 1], color = red,
                alpha = 0.3)
ax.scatter(np.linspace(0, nRec - 1, nRec), probas[y == 1, 1, 0], color = red,
           s = 46, label = 'Standard fractionation')
ax.scatter(np.linspace(0, nRec - 1, nRec), probas[y == 1, 1, 1], marker = 'v',
           color = red, s = 46, label = '28 x 2.5 Gy')
ax.scatter(np.linspace(0, nRec - 1, nRec), probas[y == 1, 1, 2], marker = 's',
           color = red, s = 46, label = '20 x 3 Gy')

#ax.scatter(np.linspace(0, 66, 67), probas[noRec, 1, 1], color = blue)
#ax.fill_between([0, 67], [0, 0], [1, 1], color = green, alpha = 0.3)
#ax.fill_between([67, 75], [0, 0], [1, 1], color = red, alpha = 0.3)
ax.set(ylabel = 'Probability', xticks = [])
ax.legend(loc = 'upper left')


#%%
fig, ax = plt.subplots()

ax.scatter(np.linspace(0, nNoRec - 1, nNoRec), probas[y == 0, 1, 1] - probas[y == 0, 1, 0],
           marker = 'v', color = green, s = 46, label = 'Recurrence (28 x 2.5 Gy)')
ax.scatter(np.linspace(0, nNoRec - 1, nNoRec), probas[y == 0, 1, 2] - probas[y == 0, 1, 0],
           marker = 's', color = green, s = 46, label = 'Recurrence (20 x 3 Gy)')

ax.scatter(np.linspace(nNoRec, 75, nRec), probas[y == 1, 1, 1] - probas[y == 1, 1, 0],
           marker = 'v', color = red, s = 46, label = 'Recurrence (28 x 2.5 Gy)')
ax.scatter(np.linspace(nNoRec, 75, nRec), probas[y == 1, 1, 2] - probas[y == 1, 1, 0],
           marker = 's', color = red, s = 46, label = 'Recurrence (20 x 3 Gy)')

#ax.scatter(np.linspace(0, 66, 67), probas[noRec, 1, 1], color = blue)
#ax.fill_between([0, 67], [0, 0], [1, 1], color = green, alpha = 0.3)
#ax.fill_between([67, 75], [0, 0], [1, 1], color = red, alpha = 0.3)
ax.set(ylabel = 'Diff. probability', xticks = [])
ax.legend(loc = 'lower left')

#%%
fig, ax = plt.subplots()

ax.scatter(np.linspace(0, nNoRec - 1, nNoRec),
           (probas[y == 0, 1, 1] - probas[y == 0, 1, 0]) / probas[y == 0, 1, 0],
           marker = 'v', color = green, s = 46, label = 'Recurrence (28 x 2.5 Gy)')
ax.scatter(np.linspace(0, nNoRec - 1, nNoRec),
           (probas[y == 0, 1, 2] - probas[y == 0, 1, 0]) / probas[y == 0, 1, 0],
           marker = 's', color = green, s = 46, label = 'Recurrence (20 x 3 Gy)')
ax.scatter(np.linspace(nNoRec, 75, nRec),
           (probas[y == 1, 1, 1] - probas[y == 1, 1, 0]) / probas[y == 1, 1, 0],
           marker = 'v', color = red, s = 46, label = 'Recurrence (28 x 2.5 Gy)')
ax.scatter(np.linspace(nNoRec, 75, nRec),
           (probas[y == 1, 1, 2] - probas[y == 1, 1, 0]) / probas[y == 1, 1, 0],
           marker = 's', color = red, s = 46, label = 'Recurrence (20 x 3 Gy)')

#ax.scatter(np.linspace(0, 66, 67), probas[noRec, 1, 1], color = blue)
#ax.fill_between([0, 67], [0, 0], [1, 1], color = green, alpha = 0.3)
#ax.fill_between([67, 75], [0, 0], [1, 1], color = red, alpha = 0.3)
ax.set(ylabel = 'Rel. diff. probability', xticks = [])
ax.legend(loc = 'lower left')