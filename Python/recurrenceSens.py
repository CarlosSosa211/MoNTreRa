#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from SALib.sample.morris import sample
from SALib.analyze.morris import analyze
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold, StratifiedKFold, LeaveOneOut
from sklearn.model_selection import cross_val_score
from sklearn.metrics import auc, confusion_matrix, roc_curve

#%%
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
# Define the model inputs

nPar = 21
boundaries = np.loadtxt('../InputFiles/refParIntRedRTArt.dat')
#names = ['tumArea', 'tumDens', 'edgeOrder', 'tumTime', 'alphaTumG1',
#              'alphaBetaTumG1', 'alphaTumGS', 'alphaBetaTumS', 'alphaTumG2',
#              'alphaBetaTumG2', 'alphaTumM', 'alphaBetaTumM', 'alphaTumG0',
#              'alphaBetaTumG0', 'hypNecThres', 'DO2', 'VmaxO2', 'KmO2',
#              'pO2NormVes']
#
#names = ['tumArea', 'tumDens', 'vascDens', 'edgeOrder', 'tumTime', 'alphaTumG1',
#              'alphaBetaTumG1', 'alphaTumGS', 'alphaBetaTumS', 'alphaTumG2',
#              'alphaBetaTumG2', 'alphaTumM', 'alphaBetaTumM', 'alphaTumG0',
#              'alphaBetaTumG0', 'hypNecThres', 'DO2', 'VmaxO2', 'KmO2',
#              'pO2NormVes']
names = ['tumArea', 'tumDens', 'vascDens', 'edgeOrder', 'tumTime', 'alphaTumG1',
              'alphaBetaTumG1', 'alphaTumGS', 'alphaBetaTumS', 'alphaTumG2',
              'alphaBetaTumG2', 'alphaTumM', 'alphaBetaTumM', 'alphaTumG0',
              'alphaBetaTumG0', 'hypNecThres', 'DO2', 'VmaxO2', 'KmO2',
              'pO2NormVes', 'stoc']
#names = ['edgeOrder', 'tumTime', 'alphaTumG1', 'alphaBetaTumG1', 'alphaTumGS',
#         'alphaBetaTumS', 'alphaTumG2', 'alphaBetaTumG2', 'alphaTumM',
#         'alphaBetaTumM', 'alphaTumG0', 'alphaBetaTumG0', 'hypNecThres', 'DO2',
#         'VmaxO2', 'KmO2', 'pO2NormVes', 'stoc']
problem = {
    'num_vars': nPar,
    'names': names,
    'bounds': boundaries
}

p = 20

#%%
X = sample(problem, 100, num_levels = p)

np.savetxt('../InputFiles/X.dat', X)


#%%
j = 3
X = np.loadtxt('../../Carlos/Results/Morris/100_20ParStoc_2Gy_Art/X.dat')
Y = np.loadtxt('../../Carlos/Results/Morris/100_20ParStoc_2Gy_Art/Y.res')
Si = analyze(problem, X, Y[:, j], conf_level = 0.95, print_to_console = True,
             num_levels = p)
Di = np.sqrt(Si['mu_star'] * Si['mu_star'] + Si['sigma'] * Si['sigma'])


#%%
N = 100
K = 3

clf = [LogisticRegression()]
tx = [Y[:, j]]

y = data[['bio_rec_6']].to_numpy().ravel()  

scores = np.zeros((N, K, len(tx)))
    
for j in range(N) :
        cv = StratifiedKFold(n_splits = K, shuffle = True)
        for i, x in enumerate(tx) :
            scores[j, :, i] = cross_val_score(clf[i], x, y, cv = cv,
                  scoring = 'roc_auc')
        
scoresMean = np.mean(scores, axis = 1)
scoresMean = pd.DataFrame(data = scoresMean, columns = ['im_feat.', 'ent'])        

print('Mean AUC')
print(scoresMean.mean(axis = 0))

print('Standard deviation AUC')
print(scoresMean.std(axis = 0))

print('Median AUC')
print(scoresMean.median(axis = 0))

#%%
namesL = ['$A_{tum}$', '$d_{tum}$', '$d_{end}$', '$N$', '$T_{tum}$',
          r'$\alpha_{tumG1}$', r'$\alpha/\beta_{tumG1}$', r'$\alpha_{tumS}$',
          r'$\alpha/\beta_{tumS}$', r'$\alpha_{tumG2}$',
          r'$\alpha/\beta_{tumG2}$', r'$\alpha_{tumM}$',
          r'$\alpha/\beta_{tumM}$', r'$\alpha_{tumG0}$',
          r'$\alpha/\beta_{tumG0}$', '$pO_2^{nec}$', '$D^{O_2}$',
          '$V_{max}^{O_2}$', '$K_M^{O_2}$','$pO_2^{preEnd}$', '$s$']
color = [brown, brown, brown, blue, blue, darkBlue, darkBlue, darkBlue, darkBlue,
         darkBlue, darkBlue, darkBlue, darkBlue, darkBlue, darkBlue, darkGreen,
         lightGreen, lightGreen, lightGreen, lightGreen, gray]

marker = ['o', 'v', '^', 'o', 'v', 'o', 'v', '^', '<', '>', 's', 'p', 'P', '*',
          'D', 'o', 'o', 'v', '^', '<', 'o']
#%%
ind = np.arange(nPar)
indSort = np.argsort(-Di)
DiSort = Di[indSort]
namesLSort = []
posStoc = np.where(indSort == 20)[0][0]

DiSortNoStoc = np.delete(DiSort, posStoc)
indSortNoStoc = np.delete(indSort, posStoc)
nParNoStoc = nPar - 1

fig, ax = plt.subplots() 
for i in range(nParNoStoc) :
    ax.bar(ind[i], DiSortNoStoc[i], facecolor = color[indSortNoStoc[i]],
           edgecolor = 'k')
    namesLSort[len(namesLSort):] = [namesL[indSortNoStoc[i]]]

ax.plot([posStoc - 0.5, posStoc - 0.5], [0, 1.1 * np.max(DiSort)], '--',
         color = gray, linewidth = 4)


ax.set(ylabel = 'Euclidean distance to the origin')
plt.xticks(ind[:-1], namesLSort, rotation = 90)

#%%
fig, ax = plt.subplots() 
for i in range(nParNoStoc) :
    ax.scatter(Si['mu_star'][indSortNoStoc[i]], Si['sigma'][indSortNoStoc[i]],
               s = 76, marker = marker[indSortNoStoc[i]],
               facecolor = color[indSortNoStoc[i]])
    namesLSort[len(namesLSort):] = [namesL[indSortNoStoc[i]]]

ax.set(xlabel = '$\mu^*$', ylabel = '$\sigma$', xlim = [0, None],
       ylim = [0, None])
ax.legend(namesLSort, bbox_to_anchor = (1.05, 1.0), loc = 'upper left',
          ncol = 2, handletextpad = 0.1, columnspacing = 0.3)

minMax = np.min([np.max(Si['mu_star']), np.max(Si['sigma'])])
ax.plot([0, 1.1 * minMax], [0, 1.1 * minMax], '--', color = gray,
        linewidth = 4)