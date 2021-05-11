#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
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
from scipy.special import xlogy
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold, StratifiedKFold, LeaveOneOut
from sklearn.model_selection import cross_val_score
from sklearn.metrics import auc, confusion_matrix, roc_curve
import radiomics as rad

#%%
plt.close('all')
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
white = [1., 1., 1.]

#%%
data = pd.read_csv('../../Carlos/Results/Recurrence/simp/rec_summary_8wTumVol_sph.csv')

#%%
path = '../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120/'

nTissues = 76
nRep = 13
N = 361
dx = 6
#rec = [9, 10, 18, 39, 41, 47, 62, 69, 72]
tRec = [9, 10, 18, 25, 26, 31, 37, 39, 41, 47, 50, 51, 57, 59, 62, 66, 69, 72]
nRec = len(tRec)
nNoRec = nTissues - len(tRec)
iRec = 0
iNoRec = 0
rec = np.zeros(nTissues)
rec[tRec] = 1 
tNoRec = np.where(rec == 0)[0].tolist()
tumVolK = np.zeros((N, 2, nRep))
tumVol = np.zeros((N, nTissues))

temp2 = np.zeros((N, 2))
temp4 = np.zeros((N, 4))
temp6 = np.zeros((N, 6))

tumDensK = np.zeros((N, nRep, nTissues))
vascDensK = np.zeros((N, nRep, nTissues))
deadDensK = np.zeros((N, nRep, nTissues))
G0DensK = np.zeros((N, nRep, nTissues))

k = 3
kern = np.ones(2 * k + 1) / (2 * k + 1)

for i in range(nTissues) :   
    for k in range(nRep) :
        temp2 = np.loadtxt(path + 'rep' + str(k) + '/tumDens_' +  str(i + 1) +
                           '.res')
        tumDensK[:, k, i] += temp2[:, 1]
        
        temp4 = np.loadtxt(path + 'rep' + str(k) + '/vascDens_' + str(i + 1) +
                           '.res')
        vascDensK[:, k, i] += temp4[:, 1]
        
        temp2 = np.loadtxt(path + 'rep' + str(k) + '/deadDens_' + str(i + 1) +
                           '.res')
        deadDensK[:, k, i] += temp2[:, 1]
        
        temp6 = np.loadtxt(path + 'rep' + str(k) + '/cycle_' + str(i + 1) +
                           '.res')
        G0DensK[:, k, i] += temp6[:, 5]

tumDensK /= 100. 
vascDensK /= 100.
deadDensK /= 100.
healDensK = 1. - tumDensK - vascDensK - deadDensK
G0DensK /= 100.

t = temp2[:, 0]
tW = t / (24 * 7)

healDens = np.median(healDensK, axis = 1)
tumDens = np.median(tumDensK, axis = 1)
vascDens = np.median(vascDensK, axis = 1)
deadDens = np.median(deadDensK, axis = 1)
G0Dens = np.median(G0DensK, axis = 1)

#%%
#entropy = - (xlogy(healDens, healDens) + xlogy(tumDens, tumDens) +
#             xlogy(vascDens, vascDens) + xlogy(deadDens, deadDens))

#entropy = G0Dens * tumDens
#entropy = jointAve / jointAve[0, :]
entropy = invDiff
meanEntRec = np.median(entropy[:, tRec], axis = 1)
stdEntRec = np.std(entropy[:, tRec], axis = 1)

meanEntNoRec = np.median(entropy[:, tNoRec], axis = 1)
stdEntNoRec = np.std(entropy[:, tNoRec], axis = 1)

wn = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
tn = wn * 28 #eight w

wEnt = entropy[tn, :]

meanWEntRec = meanEntRec[tn]
stdWEntRec = stdEntRec[tn]

meanWEntNoRec = meanEntNoRec[tn]
stdWEntNoRec = stdEntNoRec[tn]

wIntEnt = np.zeros((len(tn), nTissues))

for i in range(len(tn)) :
    wIntEnt[i, :] = np.trapz(entropy[:tn[i], :], dx = dx, axis = 0)

meanWIntEntRec = np.median(wIntEnt[:, tRec], axis = 1)
stdWIntEntRec = np.std(wIntEnt[:, tRec], axis = 1)

meanWIntEntNoRec = np.median(wIntEnt[:, tNoRec], axis = 1)
stdWIntEntNoRec = np.std(wIntEnt[:, tNoRec], axis = 1)

pvalueEnt = np.ones(N)
for i in range(3, N) :
    _, pvalueEnt[i] = stats.mannwhitneyu(entropy[i, tRec].ravel(),
                 entropy[i, tNoRec].ravel())
    
#%%
N = 361
fig, ax = plt.subplots() 

pltNoRec = ax.plot(tW[:N], entropy[:N, tNoRec], color = green, linewidth = 6,
                   alpha = 0.5, label = 'No biochemical recurrence')
plt.setp(pltNoRec[1:], label="_")
pltRec = ax.plot(tW[:N], entropy[:N, tRec], color = red, linewidth = 6,
                 alpha = 0.5, label = 'Biochemical recurrence')
plt.setp(pltRec[1:], label="_")
ax.set(xlabel = 't (weeks)', ylabel = 'Joint average',
       xlim = [0, 12.1])
ax.set_xticks(np.linspace(0, 12, 13))
ax.legend(loc = 'upper left')

#%%
fig, ax = plt.subplots() 
N = 361

ax.plot(tW[:N], meanEntRec[:N], color = red, linewidth = 6)
ax.fill_between(tW[:N], meanEntRec[:N] - stdEntRec[:N],
                meanEntRec[:N] + stdEntRec[:N], color = red,
                alpha = 0.1)
ax.plot(tW[:N], meanEntNoRec[:N], color = green, linewidth = 6)
ax.fill_between(tW[:N], meanEntNoRec[:N] - stdEntNoRec[:N],
                meanEntNoRec[:N] + stdEntNoRec[:N], color = green,
                alpha = 0.1)

#ax.plot(meanEntRec[:, 0], pvalue, color = gray, linewidth = 6)

upLim = 1.1 * max(meanEntRec[2:N] + stdEntRec[2:N]) * np.ones(N - 2)
lowLim = 1.1 * min(meanEntRec[2:N] - stdEntRec[2:N]) * np.ones(N - 2)
ax.fill_between(tW[2:N], lowLim, upLim, where = pvalueEnt[2:N] <= 0.001,
                color = gray, linewidth = 0, alpha = 0.5)
#ax.text(3.5, 1, 'p $\leq$ 0.001')
#
#ax.plot([8, 8], [0, 1.1], '--', color = black, linewidth = 6)
#ax.text(6.5, 1, 'Prediction 4', bbox = dict(facecolor = greenBlue + [0.5]))

ax.set(xlabel = 't (weeks)', ylabel = 'Joint average',
       xlim = [0, 12.1])
ax.set_xticks(np.linspace(0, 12, 13))
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'],
          loc = 'lower right')

#%%
fig, ax = plt.subplots() 

tcolor = [greenBlue, lightRed, green, red, lightGreen, redOrange]

wRec = np.tile(rec, len(wn)).reshape(nTissues * len(wn), 1)

twn = np.zeros((len(wn) * nTissues, 1))
for i in range(len(wn)) :
    twn[i * nTissues: (i + 1) * nTissues] = wn[i]

wEnt_ = np.reshape(wEnt, len(twn), order = 'C')

dataEnt = np.concatenate((wRec, twn, wEnt_.reshape((len(wEnt_), 1))),
                      axis = 1)
wEntDF = pd.DataFrame(data = dataEnt, columns = ['bio_rec', 'w', 'ent'])

ax = sns.boxplot(data = wEntDF, x = 'w', y = 'ent', hue = 'bio_rec',
                 orient = 'v', palette = tcolor)

#%%
i = 8
ind = np.argsort(wEnt[i, :])

fig, ax = plt.subplots()
for j, ii in enumerate(ind) :
    if ii in tNoRec :
        ax.scatter(j, wEnt[i, ii], color = green)
    else : 
        ax.scatter(j, wEnt[i, ii], color = red)
ax.grid()

#%%
plt.rcParams.update({'font.size': 32})
N = 100
K = 3

w = 8

clf = [LogisticRegression(), LogisticRegression(), LogisticRegression(),
       LogisticRegression(), LogisticRegression(), LogisticRegression()]
#       LogisticRegression()]
#clf = [RandomForestClassifier(), LogisticRegression()]
#tx = [data[['ADC_med', 'T2w_diff_var_mean', 'tum_vol',
#            'T2w_contrast_mean']],
#    pd.concat([data[['ADC_med', 'T2w_diff_var_mean', 'tum_vol',
#            'T2w_contrast_mean']],
#    wEntDF.loc[wEntDF['w'] == w][['ent']].reset_index(drop = True)], axis = 1)]
tx = [data[['ADC_med', 'T2w_diff_var_mean', 'tum_vol',
            'T2w_contrast_mean']],
      data[['ADC_med', 'T2w_diff_var_mean', 'tum_vol',
            'T2w_contrast_mean', 'sphericity']],
      data[['tum_vol', 'T2w_ave', 'ADC_ave']],
      data[['tum_vol', 'T2w_ave', 'ADC_ave', 'sphericity']],
      data[['init_tum_area']],
      data[['init_tum_area', 'sphericity']]]
#      wEntDF.loc[wEntDF['w'] == w][['ent']]]

y = data[['bio_rec_6']].to_numpy().ravel()  

scores = np.zeros((N, K, len(tx)))
    
for j in range(N) :
        cv = StratifiedKFold(n_splits = K, shuffle = True)
        for i, x in enumerate(tx) :
            scores[j, :, i] = cross_val_score(clf[i], x, y, cv = cv,
                  scoring = 'roc_auc')
        
scoresMean = np.mean(scores, axis = 1)
scoresMean = pd.DataFrame(data = scoresMean, columns = ['khem_feat.', 
                                                        'khem_feat_sph.',
                                                        'im_feat.',
                                                        'im_feat_sph.',
                                                        'init_tum_area',
                                                        '8w_tum_area'])
#                                                        '8w_joint_ave'])        

print('Mean AUC')
print(scoresMean.mean(axis = 0))

print('Standard deviation AUC')
print(scoresMean.std(axis = 0))

print('Median AUC')
print(scoresMean.median(axis = 0))

#%%
plt.close('all')

tcolor = [red, orange, blue, green, greenBlue, darkPurple, redOrange, redPurple]
ax = sns.boxplot(data = scoresMean, orient = 'v', palette = tcolor)

add_stat_annotation(ax, data = scoresMean,
                    box_pairs = [('khem_feat.', 'khem_feat_sph.'),
                                 ('im_feat.', 'im_feat_sph.')],
                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
                    line_offset = -0.05, verbose = 2)

ax.set(ylabel = 'AUC', ylim = [0.6, 1])
ax.set_yticks([0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_xticklabels([], ha = 'center')

#%%
plt.rcParams.update({'font.size': 32})  
N = 100
K = 3

w0 = 0
w1 = 8

#clf = [RandomForestClassifier(), LogisticRegression(), LogisticRegression()]
#
#tx = [data[['ADC_med', 'max_tum_area', 'T2w_diff_var_mean', 'tum_vol',
#            'T2w_contrast_mean']],
#    data[['tum_vol', 'ADC_ave', 'T2w_ave']],
#    wEntDF.loc[wEntDF['w'] == w][['ent']]]

clf = [RandomForestClassifier(), RandomForestClassifier(),
       LogisticRegression(), RandomForestClassifier(),
       LogisticRegression(), LogisticRegression(), LogisticRegression()]

tx = [data[['ADC_med', 'max_tum_area', 'T2w_diff_var_mean', 'tum_vol',
            'T2w_contrast_mean']],
    data[['tum_vol', 'ADC_ave', 'T2w_ave']],
    data[['tum_vol', 'ADC_ave', 'T2w_ave']],
    data[['tum_area_from_vol', 'dens_ADCT2w']],
    data[['tum_area_from_vol', 'dens_ADCT2w']],
    wEntDF.loc[wEntDF['w'] == w0][['ent']],
    wEntDF.loc[wEntDF['w'] == w1][['ent']]]
y = wEntDF.loc[wEntDF['w'] == w0]['bio_rec'].to_numpy().ravel()  


mean_fpr = np.linspace(0, 1, 500)
mean_tprK = np.zeros((N, 500, len(tx)))
auc_roc = np.zeros((N, len(tx))) 

intercept = np.zeros((N, len(tx)))
coef = []
mean_coef = []

for i, x in enumerate(tx):
    coef[len(coef):] = [np.zeros((N, len(x.columns)))]
    mean_coef[len(mean_coef):] = [np.zeros(len(x.columns))]

for j in range(N) :
    cv = StratifiedKFold(n_splits = K, shuffle = True)    
    for i, x in enumerate(tx) :
        x = x.to_numpy()
        for k, (train, test) in enumerate(cv.split(x, y)) :
#            xS, yS = SMOTEENN().fit_sample(x[train], y[train])
            xS, yS = x[train], y[train]
            xtest, ytest = x[test], y[test]
            clf[i].fit(xS, yS)
            probas = clf[i].predict_proba(xtest)
            fpr, tpr, _ = roc_curve(ytest, probas[:, 1]) 
            mean_tprK[j, :, i] += np.interp(mean_fpr, fpr, tpr)
            mean_tprK[j, 0, i] = 0.0
#            intercept[j, i] += clf[i].intercept_
#            coef[i][j, :] += clf[i].coef_.ravel()
            

intercept = intercept / K
mean_intercept = np.mean(intercept, axis = 0)

for i in range (len(tx)):
    coef[i] /= K 
    mean_coef[i] = np.mean(coef[i], axis = 0)
       
mean_tprK /= K
mean_tprK[:, -1, :] = 1.0
auc_roc = np.trapz(mean_tprK, mean_fpr, axis = 1) 
    
mean_tprN = np.mean(mean_tprK, axis = 0)
std_tprN = np.std(mean_tprK, axis = 0, ddof = 1)

level = 0.95
dof = K - 1

fig, ax_roc = plt.subplots() 
#tcolor = [blue, 'tab:gray', brown]

tcolor = [darkPurple, darkBlue, orangePurple, orange, greenBlue, darkPurple,
          redOrange, redPurple]

for i in range(len(tx)) :
    ax_roc.plot(mean_fpr, mean_tprN[:, i], linewidth = 6, color = tcolor[i])
    ci = stats.t.interval(level, dof, mean_tprN[:, i], std_tprN[:, i] / np.sqrt(N))
    ax_roc.fill_between(mean_fpr, ci[0], ci[1], color = tcolor[i],
                        alpha = 0.1)
#color = tcolor[i]
ax_roc.plot([0, 1], [0, 1], '--', color = 'tab:gray', linewidth = 6)
ax_roc.set(xlabel = 'FPR', ylabel = 'TPR', ylim = [0, 1])
mean_auc_roc = np.mean(auc_roc, axis = 0)
std_auc_roc = np.std(auc_roc, axis = 0)

    
#ax_roc.legend(["Khemara parameters\n(AUC = %0.2f $\pm$ %0.2f)" %
#           (mean_auc_roc[0], std_auc_roc[0]),
#           "Imaging parameters\n(AUC = %0.2f $\pm$ %0.2f)" %
#           (mean_auc_roc[1], std_auc_roc[1]),
#           "Entropy\n(AUC = %0.2f $\pm$ %0.2f)" %
#           (mean_auc_roc[2], std_auc_roc[2])], loc = 'lower right')
    
ax_roc.legend(["Khem features (RF) \n(AUC = %0.2f $\pm$ %0.2f)" %
           (mean_auc_roc[0], std_auc_roc[0]),
           "Imaging parameters (RF) \n(AUC = %0.2f $\pm$ %0.2f)" %
           (mean_auc_roc[1], std_auc_roc[1]),
           "Imaging parameters (LR) \n(AUC = %0.2f $\pm$ %0.2f)" %
           (mean_auc_roc[2], std_auc_roc[2]),
           "Tissue parameters (RF) \n(AUC = %0.2f $\pm$ %0.2f)" %
           (mean_auc_roc[3], std_auc_roc[3]),
           "Tissue parameters (LR) \n(AUC = %0.2f $\pm$ %0.2f)" %
           (mean_auc_roc[4], std_auc_roc[4]),
           "Tumor volume 0 weeks (LR) \n(AUC = %0.2f $\pm$ %0.2f)" %
           (mean_auc_roc[5], std_auc_roc[5]),
           "Tumor volume 8 weeks (LR) \n(AUC = %0.2f $\pm$ %0.2f)" %
           (mean_auc_roc[6], std_auc_roc[6])], loc = 'lower right')

#%%
N = 225

fig, ax = plt.subplots() 
for i in tNoRec:
    ax.plot(tW[:N], entropy[:N, i], color = green, linewidth = 6, alpha = 0.5)
    
for i in tRec:
    ax.plot(tW[:N], entropy[:N, i], color = red, linewidth = 6, alpha = 0.5)
    
#%%
path = '../../Carlos/Results/Recurrence/simp/'
nTissues = 76
nRep = 1
N = 361
dx = 6
tRec = [9, 10, 18, 25, 26, 31, 37, 39, 41, 47, 50, 51, 57, 59, 62, 66, 69, 72]
nRec = len(tRec)
nNoRec = nTissues - len(tRec)
iRec = 0
iNoRec = 0
rec = np.zeros(nTissues)
rec[tRec] = 1 
tNoRec = np.where(rec == 0)[0].tolist()

tissueDim = np.zeros((nTissues, 2), dtype = int)

state = []
pO2 = []

for i in range(nTissues) :
    temp = np.loadtxt(path + '/tissueDim/tissueDimRecSphe_' + str(i + 1) + '.dat')
    tissueDim[i, :] = temp[:2]
    state[len(state):] = [np.zeros([tissueDim[i, 0], tissueDim[i, 1], N, nRep],
          dtype = int)]
    pO2[len(pO2):] = [np.zeros([tissueDim[i, 0], tissueDim[i, 1], N, nRep],
        dtype = float)]

for i in range(nTissues) :   
    for k in range(0, 1) :
        temp = np.loadtxt(path + 'TTum330_alphaG1120_sph/rep' + str(k) +
                          '/state_' + str(i + 1) + '.res')
        state[i][:, :, :, k] = np.reshape(np.transpose(temp),
             (tissueDim[i, 0], tissueDim[i, 1], N))
        
        temp = np.loadtxt(path + 'TTum330_alphaG1120_sph/rep' + str(k) +
                          '/pO2_' + str(i + 1) + '.res')
        pO2[i][:, :, :, k] = np.reshape(np.transpose(temp),
           (tissueDim[i, 0], tissueDim[i, 1], N))
        
t = np.arange(0, N * 6, 6)
tW = t / (24 * 7)
        
#%%
statePO2 = []
sumStatePO2K = np.zeros((N, nRep, nTissues))
maxPO2 = 42.

for i in range(nTissues) :
    statePO2[len(statePO2):] = [np.zeros([tissueDim[i, 0], tissueDim[i, 1], N,
             nRep], dtype = float)]
    
for i in range(nTissues) :   
    for k in range(0, 1) :        
        statePO2[i][:, :, :, k] = np.logical_or(state[i][:, :, :, k] == 2,
                state[i][:, :, :, k] == 3) * (1. - pO2[i][:, :, :, k] / maxPO2)
    
        for j in range(N) :
            sumStatePO2K[j, k, i] = np.sum(statePO2[i][:, :, j, k])

sumStatePO2 = np.mean(sumStatePO2K, axis = 1)
sumStatePO2 /= sumStatePO2[1, :]

#%%
path = '../../Carlos/Figures/Recurrence/TissueMaps/TTum330_alphaG1120_sph/'
cmap = ListedColormap([white, blue, darkBlue, darkRed, brown])
for i in range(nTissues) :  
    plt.imsave(path + '8w_tissue' + str(i + 1) + '.png',
               state[i][:, :, 224, 0], vmin = 1, vmax = 5, cmap = cmap)
    
    
#%%
path = '../../Carlos/Figures/Recurrence/TissueMaps/TTum330_alphaG1120/'
cmap = ListedColormap([white, blue, darkBlue, darkRed, brown])
for i in range(nTissues) :  
    plt.imsave(path + '8w_tissue' + str(i + 1) + '.png',
               state[i][:, :, 0, 0], vmin = 1, vmax = 5, cmap = cmap)
    
#%%
path = '../../Carlos/Results/Recurrence/simp/'
nTissues = 76
nRep = 5
N = 361
dx = 6
tRec = [9, 10, 18, 25, 26, 31, 37, 39, 41, 47, 50, 51, 57, 59, 62, 66, 69, 72]
nRec = len(tRec)
nNoRec = nTissues - len(tRec)
iRec = 0
iNoRec = 0
rec = np.zeros(nTissues)
rec[tRec] = 1 
tNoRec = np.where(rec == 0)[0].tolist()

tissueDim = np.zeros((nTissues, 2), dtype = int)

state = []

for i in range(nTissues) :
    temp = np.loadtxt(path + '/tissueDim/tissueDimRec_' + str(i + 1) + '.dat')
    tissueDim[i, :] = temp[:2]
    state[len(state):] = [np.zeros([tissueDim[i, 0], tissueDim[i, 1], N, nRep],
          dtype = int)]

for i in range(nTissues) :   
    for k in range(10, 15) :
        temp = np.loadtxt(path + 'TTum330_alphaG1120/rep' + str(k) +
                          '/state_' + str(i + 1) + '.res')
        state[i][:, :, :, k - 10] = np.reshape(np.transpose(temp),
             (tissueDim[i, 0], tissueDim[i, 1], N))
        
t = np.arange(0, N * 6, 6)
tW = t / (24 * 7)

#%%       
windowSize = 30
windowSize2 = windowSize // 2
sWindow = windowSize * windowSize
window = np.zeros((windowSize, windowSize, N))
healDensK = np.zeros((N, nRep, nTissues))
tumDensK = np.zeros((N, nRep, nTissues))
vascDensK = np.zeros((N, nRep, nTissues))
deadDensK = np.zeros((N, nRep, nTissues))

for i in range(nTissues) :   
    nRowLow = tissueDim[i, 0] // 2 - windowSize2
    nRowUp = tissueDim[i, 0] // 2 + windowSize2
    nColLow = tissueDim[i, 1] // 2 - windowSize2
    nColUp = tissueDim[i, 1] // 2 + windowSize2
    
    for k in range(10, 15) :      
        window = state[i][nRowLow:nRowUp, nColLow:nColUp, :, k - 10]
        healDensK[:, k - 10, i] = np.count_nonzero(window == 1,
                 axis = (0, 1)) / sWindow
        tumDensK[:, k - 10, i] = np.count_nonzero(np.logical_or(window == 2,
                window == 3), axis = (0, 1))/ sWindow
#        tumDensK[:, k - 10, i] = np.count_nonzero(window == 2,
#                axis = (0, 1))/ sWindow
        vascDensK[:, k - 10, i] =  np.count_nonzero(np.logical_or(window == 4,
                window == 5), axis = (0, 1))/ sWindow
        deadDensK[:, k - 10, i] = np.count_nonzero(np.logical_or(window == 6,
                window == 7), axis = (0, 1))/ sWindow
    

healDens = np.mean(healDensK, axis = 1)
tumDens = np.mean(tumDensK, axis = 1)
vascDens = np.mean(vascDensK, axis = 1)
deadDens = np.mean(deadDensK, axis = 1)

#healDens /= healDens[0, :]
#tumDens /= tumDens[0, :]
#vascDens /= vascDens[0, :]
#deadDens /= deadDens[0, :]

t = np.arange(0, N * 6, 6)
tW = t / (24 * 7)

#%%
N = 361
windowSize = 30
windowSize2 = windowSize // 2
sWindow = windowSize * windowSize
window = np.zeros((windowSize, windowSize, N))
jointAveK = np.zeros((N, nRep, nTissues))
jointEntK = np.zeros((N, nRep, nTissues))
contrastK = np.zeros((N, nRep, nTissues))
invDiffK = np.zeros((N, nRep, nTissues))

for i in range(nTissues) :   
    nRowLow = tissueDim[i, 0] // 2 - windowSize2
    nRowUp = tissueDim[i, 0] // 2 + windowSize2
    nColLow = tissueDim[i, 1] // 2 - windowSize2
    nColUp = tissueDim[i, 1] // 2 + windowSize2
    
    for k in range(0, 1) :      
        for j in range(0, 337, 28) :
#            window = state[i][nRowLow:nRowUp, nColLow:nColUp, j, k - 10] == 2
            img = (np.logical_or(state[i][:, :, j, k] == 2,
                                 state[i][:, :, j, k] == 3)).astype(int)
#            jointAveK[j, k, i] = jointAveAllDir(img, 2)
#            jointEntK[j, k, i] = jointEntAllDir(img, 2)
#            contrastK[j, k, i] = contrastAllDir(img, 2)
            invDiffK[j, k, i] = invDiffAllDir(img, 2)
    
#jointAve = np.mean(jointAveK, axis = 1)
#jointEnt = np.mean(jointEntK, axis = 1)
#contrast = np.mean(contrastK, axis = 1)
invDiff = np.mean(invDiffK, axis = 1)

t = np.arange(0, N * 6, 6)
tW = t / (24 * 7)

#%%
np.savetxt(path + 'TTum330_alphaG1120/jointAve.res', jointAve)
np.savetxt(path + 'TTum330_alphaG1120/jointEnt.res', jointEnt)

#%%        
def coocurrence(img, levels, direction) :
    m = np.zeros((levels, levels))
    nrow, ncol = np.shape(img)
    
    if direction == 0 :
        for i in range(nrow) :
            for j in range(ncol - 1) :
                m[img[i, j], img[i, j + 1]] += 1
                
    elif direction == 45 : 
        for i in range(1, nrow) :
            for j in range(ncol - 1) :
                m[img[i, j], img[i - 1, j + 1]] += 1
            
    elif direction == 90 :
        for i in range(1, nrow) :
            for j in range(ncol) :
                m[img[i, j], img[i - 1, j]] += 1
        
    elif direction == 135 :
        for i in range(1, nrow) :
            for j in range(1, ncol) :
                m[img[i, j], img[i - 1, j - 1]] += 1
            
    return m + np.transpose(m)

#%%        
def probCoocurrence(m) :
    return m / np.sum(m)

#%%        
def jointAve(img, levels, direction) :
    m = coocurrence(img, levels, direction) 
    p = probCoocurrence(m)
    ave = 0.0
    
#    for j in range(1, levels) : 
#            ave += m[0, j] * p [0, j]
#            
#    for i in range(1, levels) : 
#            ave += m[i, 0] * p [i, 0]
#    
#    for i in range(1, levels) :
#        for j in range(1, levels) : 
#            ave += m[i, j] * p [i, j]
            
    for i in range(levels) :
        for j in range(levels) : 
            ave += m[i, j] * p [i, j]
            
    return ave
    
#%%        
def jointEnt(img, levels, direction) :
    m = coocurrence(img, levels, direction) 
    p = probCoocurrence(m)
    ent = 0.0
    
#    for j in range(1, levels) : 
#        ent -= xlogy(m[0, j], p [0, j])
#            
#    for i in range(1, levels) :
#        ent -= xlogy(m[i, 0], p [i, 0])
#    
#    for i in range(1, levels) :
#        for j in range(1, levels) : 
#            ent -= xlogy(m[i, j], p [i, j])
#            
    for i in range(levels) :
        for j in range(levels) : 
            ent -= xlogy(m[i, j], p [i, j])
            
    return ent

#%%        
def contrast(img, levels, direction) :
    m = coocurrence(img, levels, direction) 
    p = probCoocurrence(m)
    cont = 0.0
    
#    for j in range(1, levels) : 
#        ent -= xlogy(m[0, j], p [0, j])
#            
#    for i in range(1, levels) :
#        ent -= xlogy(m[i, 0], p [i, 0])
#    
#    for i in range(1, levels) :
#        for j in range(1, levels) : 
#            ent -= xlogy(m[i, j], p [i, j])
#            
    for i in range(levels) :
        for j in range(levels) : 
            cont += (i - j) * (i - j) * p [i, j]
            
    return cont

#%%
def invDiff(img, levels, direction) :
    m = coocurrence(img, levels, direction) 
    p = probCoocurrence(m)
    idiff = 0.0
    
#    for j in range(1, levels) : 
#        ent -= xlogy(m[0, j], p [0, j])
#            
#    for i in range(1, levels) :
#        ent -= xlogy(m[i, 0], p [i, 0])
#    
#    for i in range(1, levels) :
#        for j in range(1, levels) : 
#            ent -= xlogy(m[i, j], p [i, j])
#            
    for i in range(levels) :
        for j in range(levels) : 
            idiff += p [i, j] / (1. + abs(i - j))
            
    return idiff

#%%
def jointAveAllDir(img, levels) :
    return 0.25 * (jointAve(img, levels, 0) + jointAve(img, levels, 45) +
                   jointAve(img, levels, 90) + jointAve(img, levels, 135))
    
#%%
def jointEntAllDir(img, levels) :
    return 0.25 * (jointEnt(img, levels, 0) + jointEnt(img, levels, 45) +
                   jointEnt(img, levels, 90) + jointEnt(img, levels, 135))
    

#%%
def contrastAllDir(img, levels) :
    return 0.25 * (contrast(img, levels, 0) + contrast(img, levels, 45) +
                   contrast(img, levels, 90) + contrast(img, levels, 135))

#%%    
def invDiffAllDir(img, levels) :
    return 0.25 * (invDiff(img, levels, 0) + invDiff(img, levels, 45) +
                   invDiff(img, levels, 90) + invDiff(img, levels, 135))

#%%
#entropy = - (xlogy(healDens, healDens) + xlogy(tumDens, tumDens) +
#             xlogy(vascDens, vascDens) + xlogy(deadDens, deadDens))

#entropy = healDens / tumDens
entropy = sumStatePO2

meanEntRec = np.mean(entropy[:, tRec], axis = 1)
stdEntRec = np.std(entropy[:, tRec], axis = 1)

meanEntNoRec = np.mean(entropy[:, tNoRec], axis = 1)
stdEntNoRec = np.std(entropy[:, tNoRec], axis = 1)

wn = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
tn = wn * 28 #24/6 * 7

wEnt = entropy[tn, :]

meanWEntRec = meanEntRec[tn]
stdWEntRec = stdEntRec[tn]

meanWEntNoRec = meanEntNoRec[tn]
stdWEntNoRec = stdEntNoRec[tn]

wIntEnt = np.zeros((len(tn), nTissues))

for i in range(len(tn)) :
    wIntEnt[i, :] = np.trapz(entropy[:tn[i], :], dx = dx, axis = 0)

meanWIntEntRec = np.mean(wIntEnt[:, tRec], axis = 1)
stdWIntEntRec = np.std(wIntEnt[:, tRec], axis = 1)

meanWIntEntNoRec = np.mean(wIntEnt[:, tNoRec], axis = 1)
stdWIntEntNoRec = np.std(wIntEnt[:, tNoRec], axis = 1)

pvalueEnt = np.ones(N)
for i in range(3, N) :
    _, pvalueEnt[i] = stats.mannwhitneyu(entropy[i, tRec].ravel(),
                 entropy[i, tNoRec].ravel())

#%%
fig, ax = plt.subplots() 
N = 361

ax.plot(tW[2:N], meanEntRec[2:N], color = red, linewidth = 6)
ax.fill_between(tW[2:N], meanEntRec[2:N] - stdEntRec[2:N],
                meanEntRec[2:N] + stdEntRec[2:N], color = red,
                alpha = 0.1)
ax.plot(tW[2:N], meanEntNoRec[2:N], color = green, linewidth = 6)
ax.fill_between(tW[2:N], meanEntNoRec[2:N] - stdEntNoRec[2:N],
                meanEntNoRec[2:N] + stdEntNoRec[2:N], color = green,
                alpha = 0.1)


#ax.plot(meanEntRec[:, 0], pvalue, color = gray, linewidth = 6)


upLim = 1.1 * max(meanEntRec[2:N] + stdEntRec[2:N]) * np.ones(N - 2)
lowLim = 1.1 * min(meanEntRec[2:N] - stdEntRec[2:N]) * np.ones(N - 2)
ax.fill_between(tW[2:N], lowLim, upLim, where = pvalueEnt[2:N] <= 0.001,
                color = gray, linewidth = 0, alpha = 0.5)
#ax.text(3.5, 1, 'p $\leq$ 0.001')
#
#ax.plot([8, 8], [0, 1.1], '--', color = black, linewidth = 6)
#ax.text(6.5, 1, 'Prediction 4', bbox = dict(facecolor = greenBlue + [0.5]))

ax.set(xlabel = 't (weeks)', ylabel = 'Normalized $pO_2$ n° of tumor cells',
       xlim = [0, 12.1])
ax.set_xticks(np.linspace(0, 12, 13))
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'],
          loc = 'upper right')