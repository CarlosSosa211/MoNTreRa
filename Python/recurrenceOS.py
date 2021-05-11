#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as stats
import seaborn as sns
from statannot import add_stat_annotation
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import LinearSVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import auc, confusion_matrix, roc_curve
from sklearn.model_selection import StratifiedKFold
from imblearn.over_sampling import SMOTE, SVMSMOTE, BorderlineSMOTE, KMeansSMOTE, ADASYN
from imblearn.under_sampling import RandomUnderSampler
from imblearn.combine import SMOTEENN
from imblearn.pipeline import Pipeline

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
dataOS = pd.read_csv('../../Carlos/Results/Recurrence/simp/rec_summary_8wTumVolOS2.csv')

N = 100
K = 3
x = dataOS[['tum_vol', 'ADC_ave', 'T2w_ave']].to_numpy()
dataOS.insert(loc = 0, column = 'rep', value = -1)
dataOS.insert(loc = 1, column = 'fold', value = -1)
dataOS.insert(loc = 2, column = 'test', value = -1)

y = dataOS[['bio_rec']].to_numpy().ravel() 

for j in range(N) :
    cv = StratifiedKFold(n_splits = K, shuffle = True)
    os = SMOTE()
    for k, (train_, test_) in enumerate(cv.split(x, y)) :
        xS, yS = os.fit_resample(x[train_], y[train_])
        numS = len(xS)
        rep = j * np.ones([numS, 1])
        fold = k * np.ones([numS, 1])
        test = np.zeros([numS, 1])
        syn = np.zeros([numS, 1])
        syn[len(x[train_]):] = 1
        yS = yS.reshape(numS, 1)
        datak = pd.DataFrame(np.concatenate((rep, fold, test, syn, yS, xS),
                                            axis = 1),
                                            columns = ['rep', 'fold', 'test',
                                                       'syn', 'bio_rec',
                                                       'tum_vol', 'ADC_ave',
                                                       'T2w_ave'])
        dataOS = dataOS.append(datak, ignore_index = True)
        
        numTest = len(test_)
        rep = j * np.ones([numTest, 1])
        fold = k * np.ones([numTest, 1])
        test = np.ones([numTest, 1])
        syn = np.zeros([numTest, 1])
        datak = pd.DataFrame(np.concatenate((rep, fold, test, syn, y[test_].reshape(numTest, 1),
                                             x[test_]), axis = 1),
                                            columns = ['rep', 'fold', 'test',
                                                       'syn', 'bio_rec',
                                                       'tum_vol', 'ADC_ave',
                                                       'T2w_ave'])
        dataOS = dataOS.append(datak, ignore_index = True)
        
        
ADC_ave = dataOS[['ADC_ave']]. to_numpy()
T2w_ave = dataOS[['T2w_ave']]. to_numpy()
dens_ADCnorm = - 0.56 * (ADC_ave - np.mean(ADC_ave)) / np.std(ADC_ave) 
dens_T2wnorm = - 0.39 * (T2w_ave - np.mean(T2w_ave)) / np.std(T2w_ave)

dens_ADCT2w = 0.375 * ((dens_ADCnorm - np.min(dens_ADCnorm)) / 
                      (np.max(dens_ADCnorm) - np.min(dens_ADCnorm)) +
                       (dens_T2wnorm - np.min(dens_T2wnorm)) /
                       (np.max(dens_T2wnorm) - np.min(dens_T2wnorm))) + 0.25
dataOS['dens_ADCT2w'] = dens_ADCT2w

tum_vol = dataOS[['tum_vol']]. to_numpy()
tum_area_from_vol = np.pi ** (1/3) * (0.75 * tum_vol) ** (2/3)
dataOS['tum_area_from_vol'] = tum_area_from_vol

init_tum_area = tum_area_from_vol * dens_ADCT2w
dataOS['init_tum_area'] = init_tum_area

dataOS.to_csv('../../Carlos/Results/Recurrence/simp/dataOS.csv', index = False)
dataOS.loc[dataOS['syn'] == 1][['tum_area_from_vol', 'dens_ADCT2w']].to_csv('../../Carlos/Results/Recurrence/simp/dataOSSyn.csv',
           header = False, index = False)

#%%
path = '../../Carlos/Results/Recurrence/simp/'
dataOS = pd.read_csv(path + 'dataSMOTE.csv')
dataOS.insert(loc = 0, column = 'ID', value = 0)
dataOS.insert(loc = len(dataOS.columns), column = 'init_tum_area', value = 0)
dataOS.insert(loc = len(dataOS.columns), column = '8w_tum_area', value = 0)
dataOS.insert(loc = len(dataOS.columns), column = '8w_tum_area_norm', value = 0)
dataOS.insert(loc = len(dataOS.columns), column = '8w_int_tum_area', value = 0)
dataOS.insert(loc = len(dataOS.columns), column = '8w_int_tum_area_norm', value = 0)

initTumVol = np.loadtxt(path + 'TTum330_alphaG1120_syn/initTumVol.res')
eightwTumVol = np.loadtxt(path + 'TTum330_alphaG1120_syn/8wTumVol.res')
eightwIntTumVol = np.loadtxt(path + 'TTum330_alphaG1120_syn/8wIntTumVol.res')
for i in range(76) :
    tissuei = dataOS.index[(dataOS['tum_vol'] == dataOS.at[i, 'tum_vol']) &
                        (dataOS['ADC_ave'] == dataOS.at[i, 'ADC_ave']) &
                        (dataOS['T2w_ave'] == dataOS.at[i, 'T2w_ave'])].to_list()
    dataOS.at[tissuei, 'ID'] = i + 1
    dataOS.at[tissuei, 'init_tum_area'] = initTumVol[i]
    dataOS.at[tissuei, '8w_tum_area'] = eightwTumVol[i]
    dataOS.at[tissuei, '8w_int_tum_area'] = eightwIntTumVol[i]

initTumVolSyn = np.loadtxt(path + 'TTum330_alphaG1120_syn/initTumVolSyn.res')
eightwTumVolSyn = np.loadtxt(path + 'TTum330_alphaG1120_syn/8wTumVolSyn.res')
eightwIntTumVolSyn = np.loadtxt(path + 'TTum330_alphaG1120_syn/8wIntTumVolSyn.res')

dataOS.at[dataOS.index[(dataOS['syn'] == 1)].to_list(), 'init_tum_area'] = initTumVolSyn
dataOS.at[dataOS.index[(dataOS['syn'] == 1)].to_list(), '8w_tum_area'] = eightwTumVolSyn
dataOS.at[dataOS.index[(dataOS['syn'] == 1)].to_list(), '8w_int_tum_area'] = eightwIntTumVolSyn

dataOS['8w_tum_area_norm'] = dataOS['8w_tum_area'] / dataOS['init_tum_area'] 
dataOS['8w_int_tum_area_norm'] =  dataOS['8w_int_tum_area'] / (dataOS['init_tum_area'] * 2160) 

#%%
N = 100
K = 3
logReg = LogisticRegression()
#logReg = RandomForestClassifier()
#logReg = MLPClassifier()
#logReg = GaussianProcessClassifier()
#logReg = AdaBoostClassifier()
#logReg = KNeighborsClassifier(n_neighbors = 10)
#logReg = LinearSVC()
mean_fpr = np.linspace(0, 1, 500)
mean_tprK = np.zeros((N, 500, 4))
auc_ = np.zeros((N, 4)) 

for j in range(N) :
    for k in range(K) :
        xtrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 0) & (dataOS['syn'] == 0)][['tum_vol', 'ADC_ave', 'T2w_ave']].to_numpy()
        ytrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 0) & (dataOS['syn'] == 0)]['bio_rec'].to_numpy().ravel()
        xtest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 1)][['tum_vol', 'ADC_ave', 'T2w_ave']].to_numpy()
        ytest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 1)]['bio_rec'].to_numpy().ravel()
        logReg.fit(xtrain, ytrain)
        probas = logReg.predict_proba(xtest)
        fpr, tpr, _ = roc_curve(ytest, probas[:, 1]) 
        mean_tprK[j, :, 0] += np.interp(mean_fpr, fpr, tpr)
        mean_tprK[j, 0, 0] = 0.0
        
        xtrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 0)][['tum_vol', 'ADC_ave', 'T2w_ave']].to_numpy()
        ytrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 0)]['bio_rec'].to_numpy().ravel()
        xtest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 1)][['tum_vol', 'ADC_ave', 'T2w_ave']].to_numpy()
        ytest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 1)]['bio_rec'].to_numpy().ravel()
        logReg.fit(xtrain, ytrain)
        probas = logReg.predict_proba(xtest)
        fpr, tpr, _ = roc_curve(ytest, probas[:, 1]) 
        mean_tprK[j, :, 1] += np.interp(mean_fpr, fpr, tpr)
        mean_tprK[j, 0, 1] = 0.0
        
#        xtrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
#                            (dataOS['test'] == 0)][['tum_area_from_vol', 'dens_ADCT2w']].to_numpy()
#        ytrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
#                            (dataOS['test'] == 0)]['bio_rec'].to_numpy().ravel()
#        xtest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
#                            (dataOS['test'] == 1)][['tum_area_from_vol', 'dens_ADCT2w']].to_numpy()
#        ytest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
#                            (dataOS['test'] == 1)]['bio_rec'].to_numpy().ravel()
#        logReg.fit(xtrain, ytrain)
#        probas = logReg.predict_proba(xtest)
#        fpr, tpr, _ = roc_curve(ytest, probas[:, 1]) 
#        mean_tprK[j, :, 2] += np.interp(mean_fpr, fpr, tpr)
#        mean_tprK[j, 0, 2] = 0.0
        
#        xtrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
#                            (dataOS['test'] == 0) & (dataOS['syn'] == 0)][['tum_area_from_vol', 'dens_ADCT2w']].to_numpy()
#        ytrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
#                            (dataOS['test'] == 0) & (dataOS['syn'] == 0)]['bio_rec'].to_numpy().ravel()
#        xtest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
#                            (dataOS['test'] == 1)][['tum_area_from_vol', 'dens_ADCT2w']].to_numpy()
#        ytest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
#                            (dataOS['test'] == 1)]['bio_rec'].to_numpy().ravel()
#        logReg.fit(xtrain, ytrain)
#        probas = logReg.predict_proba(xtest)
#        fpr, tpr, _ = roc_curve(ytest, probas[:, 1]) 
#        mean_tprK[j, :, 2] += np.interp(mean_fpr, fpr, tpr)
#        mean_tprK[j, 0, 2] = 0.0
#        
#        xtrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
#                            (dataOS['test'] == 0) & (dataOS['syn'] == 0)][['init_tum_area']].to_numpy()
#        ytrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
#                            (dataOS['test'] == 0) & (dataOS['syn'] == 0)]['bio_rec'].to_numpy().ravel()
#        xtest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
#                            (dataOS['test'] == 1)][['init_tum_area']].to_numpy()
#        ytest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
#                            (dataOS['test'] == 1)]['bio_rec'].to_numpy().ravel()
#        logReg.fit(xtrain, ytrain)
#        probas = logReg.predict_proba(xtest)
#        fpr, tpr, _ = roc_curve(ytest, probas[:, 1]) 
#        mean_tprK[j, :, 3] += np.interp(mean_fpr, fpr, tpr)
#        mean_tprK[j, 0, 3] = 0.0
        
        xtrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 0) & (dataOS['syn'] == 0)][['8w_tum_area_norm']].to_numpy()
        ytrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 0) & (dataOS['syn'] == 0)]['bio_rec'].to_numpy().ravel()
        xtest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 1)][['8w_tum_area_norm']].to_numpy()
        ytest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 1)]['bio_rec'].to_numpy().ravel()
        logReg.fit(xtrain, ytrain)
        probas = logReg.predict_proba(xtest)
        fpr, tpr, _ = roc_curve(ytest, probas[:, 1]) 
        mean_tprK[j, :, 2] += np.interp(mean_fpr, fpr, tpr)
        mean_tprK[j, 0, 2] = 0.0
        
        xtrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 0)][['8w_tum_area_norm']].to_numpy()
        ytrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 0)]['bio_rec'].to_numpy().ravel()
        xtest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 1)][['8w_tum_area_norm']].to_numpy()
        ytest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 1)]['bio_rec'].to_numpy().ravel()
        logReg.fit(xtrain, ytrain)
        probas = logReg.predict_proba(xtest)
        fpr, tpr, _ = roc_curve(ytest, probas[:, 1]) 
        mean_tprK[j, :, 3] += np.interp(mean_fpr, fpr, tpr)
        mean_tprK[j, 0, 3] = 0.0
            
mean_tprK /= K
mean_tprK[:, -1, :] = 1.0
auc_ = np.trapz(mean_tprK, mean_fpr, axis = 1) 
    
mean_tprN = np.mean(mean_tprK, axis = 0)
std_tprN = np.std(mean_tprK, axis = 0, ddof = 1)
level = 0.95
dof = K - 1

#%%
fig, ax = plt.subplots() 
tcolor = [red, orange, blue, green]

for i in range(4) :
    ax.plot(mean_fpr, mean_tprN[:, i], linewidth = 6, color = tcolor[i])
    ci = stats.t.interval(level, dof, mean_tprN[:, i], std_tprN[:, i] / np.sqrt(N))
    ax.fill_between(mean_fpr, ci[0], ci[1], color = tcolor[i],
                    alpha = 0.1)

ax.plot([0, 1], [0, 1], '--', color = 'tab:gray', linewidth = 6)
ax.set(xlabel = 'FPR', ylabel = 'TPR', ylim = [0, 1])
mean_auc = np.mean(auc_, axis = 0)
std_auc = np.std(auc_, axis = 0)

ax.legend(["Imaging parameters\n(AUC = %0.2f $\pm$ %0.2f)" %
           (mean_auc[0], std_auc[0]),
           "Imaging parameters\nwith SMOTE\n(AUC = %0.2f $\pm$ %0.2f)" %
           (mean_auc[1], std_auc[1]),
           "Norm. int. of tum. area\nup to 8 w.\n(AUC = %0.2f $\pm$ %0.2f)" % 
           (mean_auc[2], std_auc[2]),
           "Norm. int. of tum. area\nup to 8 w. with SMOTE\n(AUC = %0.2f $\pm$ %0.2f)" %
           (mean_auc[3], std_auc[3])],
           loc = 'lower right')

auc_ = pd.DataFrame(data = auc_, columns = ['im_feat.', 'im_feat._SMOTE',
                                            '8w_int_tum_area_norm', '8w_int_tum_area_norm_SMOTE'])
print('Mean AUC')
print(auc_.mean(axis = 0))

print('Standard deviation AUC')
print(auc_.std(axis = 0))

print('Median AUC')
print(auc_.median(axis = 0))

#%%
N = 100
K = 3
auc_i = np.zeros(10)
for i, x in enumerate(range(100, 1100, 100)) : 
    logReg = RandomForestClassifier(n_estimators = x)
    mean_fpr = np.linspace(0, 1, 500)
    mean_tprK = np.zeros((N, 500, 2))
    auc_ = np.zeros((N, 2)) 

    for j in range(N) :
        for k in range(K) :        
            xtrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 0) & (dataOS['syn'] == 0)][['8w_tum_area']].to_numpy()
            ytrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 0) & (dataOS['syn'] == 0)]['bio_rec'].to_numpy().ravel()
            xtest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 1)][['8w_tum_area']].to_numpy()
            ytest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 1)]['bio_rec'].to_numpy().ravel()
            logReg.fit(xtrain, ytrain)
            probas = logReg.predict_proba(xtest)
            fpr, tpr, _ = roc_curve(ytest, probas[:, 1]) 
            mean_tprK[j, :, 0] += np.interp(mean_fpr, fpr, tpr)
            mean_tprK[j, 0, 0] = 0.0
        
            xtrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 0)][['8w_tum_area']].to_numpy()
            ytrain = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 0)]['bio_rec'].to_numpy().ravel()
            xtest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 1)][['8w_tum_area']].to_numpy()
            ytest = dataOS.loc[(dataOS['rep'] == j) & (dataOS['fold'] == k) &
                            (dataOS['test'] == 1)]['bio_rec'].to_numpy().ravel()
            logReg.fit(xtrain, ytrain)
            probas = logReg.predict_proba(xtest)
            fpr, tpr, _ = roc_curve(ytest, probas[:, 1]) 
            mean_tprK[j, :, 1] += np.interp(mean_fpr, fpr, tpr)
            mean_tprK[j, 0, 1] = 0.0
            
    mean_tprK /= K
    mean_tprK[:, -1, :] = 1.0
    auc_ = np.trapz(mean_tprK, mean_fpr, axis = 1) 
    
    auc_i[i] = np.median(auc_, axis = 0)[1]

    print(x)
    print('Mean AUC')
    print(np.mean(auc_, axis = 0))

    print('Standard deviation AUC')
    print(np.std(auc_, axis = 0))

    print('Median AUC')
    print(np.median(auc_, axis = 0))
    
print(auc_i)

#%%
fig, ax = plt.subplots() 
ax = sns.boxplot(data = auc_, orient = 'v', palette = tcolor)
#ax = sns.swarmplot(data = scoresMean, color = ".25")
add_stat_annotation(ax, data = auc_, box_pairs = [('im_feat.','im_feat._SMOTE'),
                                                  ('8w_tum_area', '8w_tum_area_SMOTE')],
                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
                    line_offset = -0.05, verbose = 2)
#add_stat_annotation(ax, data = scoresMean,
#                    box_pairs = [('090',
#                                '120'),
#                                ('120',
#                                '154'),
#                                ('154',
#                                '223')],
#                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
#                    line_offset = -0.05, verbose = 2)
ax.set(ylabel = 'AUC', ylim = [0.5, 1])
ax.set_xticklabels([], ha = 'center')

#%%
fig, ax = plt.subplots()
plt.scatter(dataOS.loc[(dataOS['bio_rec'] == 0) & (dataOS['syn'] == 0)]['tum_area_from_vol'].to_numpy(),
            dataOS.loc[(dataOS['bio_rec'] == 0) & (dataOS['syn'] == 0)]['dens_ADCT2w'].to_numpy())
plt.scatter(dataOS.loc[(dataOS['bio_rec'] == 1) & (dataOS['syn'] == 0)]['tum_area_from_vol'].to_numpy(),
            dataOS.loc[(dataOS['bio_rec'] == 1) & (dataOS['syn'] == 0)]['dens_ADCT2w'].to_numpy())

ax.set(title = 'Density', xlabel = 'tum_area_from_vol', ylabel = 'dens_ADCT2w')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

#%%
fig, ax = plt.subplots()
plt.scatter(dataOS.loc[(dataOS['bio_rec'] == 0) & (dataOS['test'] != 1)]['tum_area_from_vol'].to_numpy(),
            dataOS.loc[(dataOS['bio_rec'] == 0) & (dataOS['test'] != 1)]['dens_ADCT2w'].to_numpy())
plt.scatter(dataOS.loc[(dataOS['bio_rec'] == 1) & (dataOS['test'] != 1)]['tum_area_from_vol'].to_numpy(),
            dataOS.loc[(dataOS['bio_rec'] == 1) & (dataOS['test'] != 1)]['dens_ADCT2w'].to_numpy())

ax.set(title = 'Density', xlabel = 'tum_area_from_vol', ylabel = 'dens_ADCT2w')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

#%%
fig, ax = plt.subplots()
plt.scatter(dataOS.loc[(dataOS['bio_rec'] == 0) & (dataOS['syn'] == 0)]['8w_tum_area'].to_numpy(),
            dataOS.loc[(dataOS['bio_rec'] == 0) & (dataOS['syn'] == 0)]['8w_tum_area'].to_numpy())
plt.scatter(dataOS.loc[(dataOS['bio_rec'] == 1) & (dataOS['syn'] == 0)]['8w_tum_area'].to_numpy(),
            dataOS.loc[(dataOS['bio_rec'] == 1) & (dataOS['syn'] == 0)]['8w_tum_area'].to_numpy())

ax.set(title = 'Density', xlabel = '8w_tum_area', ylabel = '8w_tum_area')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

#%%
fig, ax = plt.subplots()
plt.scatter(dataOS.loc[(dataOS['bio_rec'] == 0) & (dataOS['test'] != 1) & (dataOS['rep'] == 0)]['8w_tum_area'].to_numpy(),
            dataOS.loc[(dataOS['bio_rec'] == 0) & (dataOS['test'] != 1) & (dataOS['rep'] == 0)]['8w_tum_area'].to_numpy())
plt.scatter(dataOS.loc[(dataOS['bio_rec'] == 1) & (dataOS['test'] != 1) & (dataOS['rep'] == 0)]['8w_tum_area'].to_numpy(),
            dataOS.loc[(dataOS['bio_rec'] == 1) & (dataOS['test'] != 1) & (dataOS['rep'] == 0)]['8w_tum_area'].to_numpy())

ax.set(title = 'Density', xlabel = 'tum_area_from_vol', ylabel = 'dens_ADCT2w')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

#%%
fig, ax = plt.subplots()
sns.distplot(dataOS.loc[(dataOS['bio_rec'] == 1) & (dataOS['syn'] == 0) & (dataOS['rep'] == 0)]['8w_tum_area'].to_numpy())
sns.distplot(dataOS.loc[(dataOS['bio_rec'] == 0) & (dataOS['syn'] == 0) & (dataOS['rep'] == 0)]['8w_tum_area'].to_numpy())
ax.set(title = 'Without SMOTE', xlabel = '8w_tum_area', ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(dataOS.loc[(dataOS['bio_rec'] == 1) & (dataOS['syn'] == 0) & (dataOS['rep'] == 0)]['8w_tum_area'].to_numpy(),
                      dataOS.loc[(dataOS['bio_rec'] == 0) & (dataOS['syn'] == 0) & (dataOS['rep'] == 0)]['8w_tum_area'].to_numpy())
print('Without SMOTE', p)

#%%
fig, ax = plt.subplots()
sns.distplot(dataOS.loc[(dataOS['bio_rec'] == 1) & (dataOS['test'] != 1) & (dataOS['rep'] == 0)]['8w_tum_area'].to_numpy())
sns.distplot(dataOS.loc[(dataOS['bio_rec'] == 0) & (dataOS['test'] != 1) & (dataOS['rep'] == 0)]['8w_tum_area'].to_numpy())
ax.set(title = 'With SMOTE', xlabel = '8w_tum_area', ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(dataOS.loc[(dataOS['bio_rec'] == 1) & (dataOS['test'] != 1) & (dataOS['rep'] == 0)]['8w_tum_area'].to_numpy(),
                      dataOS.loc[(dataOS['bio_rec'] == 0) & (dataOS['test'] != 1) & (dataOS['rep'] == 0)]['8w_tum_area'].to_numpy())
print('With SMOTE', p)

#%%
fig, ax = plt.subplots()
plt.scatter(dataOS.loc[(dataOS['bio_rec'] == 0) & (dataOS['test'] != 1)]['T2w_ave'].to_numpy(),
            dataOS.loc[(dataOS['bio_rec'] == 0) & (dataOS['test'] != 1)]['ADC_ave'].to_numpy())
plt.scatter(dataOS.loc[(dataOS['bio_rec'] == 1) & (dataOS['test'] != 1)]['T2w_ave'].to_numpy(),
            dataOS.loc[(dataOS['bio_rec'] == 1) & (dataOS['test'] != 1)]['ADC_ave'].to_numpy())

ax.set(title = 'Density', xlabel = 'T2w', ylabel = 'ADC')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

#%%
fig, ax = plt.subplots()
plt.scatter(dataOS.loc[dataOS['bio_rec'] == 0]['tum_area_from_vol'].to_numpy(),
            dataOS.loc[dataOS['bio_rec'] == 0]['dens_ADCT2w'].to_numpy())
plt.scatter(dataOS.loc[dataOS['bio_rec'] == 1]['tum_area_from_vol'].to_numpy(),
            dataOS.loc[dataOS['bio_rec'] == 1]['dens_ADCT2w'].to_numpy())

ax.set(title = 'Density', xlabel = 'tum_area_from_vol', ylabel = 'dens_ADCT2w')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

#%%
data = pd.read_csv('../../Carlos/Results/Recurrence/simp/rec_summary_8wTumVol.csv')

#%%
plt.rcParams.update({'font.size': 32})
N = 1000
K = 3
logReg = LogisticRegression()
#logReg = LinearSVC()
#logReg = RandomForestClassifier()
#logReg = MLPClassifier(activation = 'logistic', max_iter = 1000)

#tx = data[['ADC_med', 'max_tum_area', 'T2w_diff_var_mean', 'tum_vol',
#            'T2w_contrast_mean']]

tx = [data[['tum_vol', 'ADC_ave', 'T2w_ave']],
      data[['tum_area_from_vol', 'dens_ADCT2w']],
      data[['init_tum_area']],
      data[['TTum330_alphaG1120_norm']]]

#tx = [data[['TTum260_alphaG1120']],
#      data[['TTum330_alphaG1120']],
#      data[['TTum360_alphaG1120']],
#      data[['TTum400_alphaG1120']]]

#tx = [data[['TTum330_alphaG1090']],
#      data[['TTum330_alphaG1120']],
#      data[['TTum330']],
#      data[['TTum330_alphaG1223']]]

y = data[['bio_rec']].to_numpy().ravel() 

pipeline = [Pipeline(steps = [('model', logReg)]),
            Pipeline(steps = [('over', SMOTEENN()), ('model', logReg)])]

scores = np.zeros((N, K, 2 * len(tx)))
    
for j in range(N) :
        cv = StratifiedKFold(n_splits = K, shuffle = True)
        for i, x in enumerate(tx) :
            scores[j, :, 2 * i] = cross_val_score(pipeline[0], x, y, cv = cv,
                  scoring = 'roc_auc')
            scores[j, :, 2 * i + 1] = cross_val_score(pipeline[1], x, y, cv = cv,
                  scoring = 'roc_auc')
        
scoresMean = np.mean(scores, axis = 1)
scoresMean = pd.DataFrame(data = scoresMean,
                          columns = ['im_feat.', 'im_feat_OS', 'tissue_feat.',
                                     'tissue_feat_OS', 'init_tum_area',
                                     'init_tum_area_OS', '8w_tum_area_norm',
                                     '8w_tum_area_norm_OS'])
#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['260', '330', '360', '400'])
#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['090', '120', '154', '223'])

print('Mean AUC')
print(scoresMean.mean(axis = 0))

print('Standard deviation AUC')
print(scoresMean.std(axis = 0))

print('Median AUC')
print(scoresMean.median(axis = 0))

#%%
#tcolor = [darkPurple, darkPurple, red, red, orange, orange, green, green]
tcolor = [darkPurple, darkBlue, orangePurple, orange, greenBlue]
fig, ax = plt.subplots() 

scoresPlot = scoresMean[['im_feat.', 'im_feat_OS', 'tissue_feat.',
                       'init_tum_area', '8w_tum_area_norm']]
ax = sns.boxplot(data = scoresPlot, orient = 'v', palette = tcolor,
                 linewidth = 2, fliersize = 6)
#ax = sns.swarmplot(data = scoresMean, color = ".25")
add_stat_annotation(ax, data = scoresPlot, box_pairs = [('im_feat.','im_feat_OS'),
                                                  ('im_feat.', 'tissue_feat.'),
                                                  ('tissue_feat.', 'init_tum_area'),
                                                  ('init_tum_area', '8w_tum_area_norm')],
                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
                    linewidth = 2, line_offset = -0.05, verbose = 2)

ax.set(ylabel = 'AUC', ylim = [0.57, 1])
ax.set_xticklabels(['Prediction 1', "Prediction 1*", 'Prediction 2',
                    'Prediction 3', 'Prediction 4'], ha = 'center')
ax.set_yticks([0.6, 0.7, 0.8, 0.9, 1.0])

#%%
plt.close('all')
plt.rcParams.update({'font.size': 32})  
N = 1000
K = 3

clf = [LogisticRegression(), LogisticRegression(), LogisticRegression(),
       LogisticRegression(), LogisticRegression(), LogisticRegression(),
       LogisticRegression(), LogisticRegression()]

tx = [data[['tum_vol', 'ADC_ave', 'T2w_ave']],
      data[['tum_vol', 'ADC_ave', 'T2w_ave']],
      data[['tum_area_from_vol', 'dens_ADCT2w']],
      data[['tum_area_from_vol', 'dens_ADCT2w']],
      data[['init_tum_area']], 
      data[['init_tum_area']],
      data[['TTum330_alphaG1120_norm']],
      data[['TTum330_alphaG1120_norm']]]

y = data[['bio_rec']].to_numpy().ravel()  

mean_fpr = np.linspace(0, 1, 500)
mean_tprK = np.zeros((N, 500, len(tx)))
auc_roc = np.zeros((N, len(tx))) 

intercept = np.zeros((N, len(tx)))
coef = []
mean_coef = []
sme = SMOTEENN(random_state = 42)

for i, x in enumerate(tx):
    coef[len(coef):] = [np.zeros((N, len(x.columns)))]
    mean_coef[len(mean_coef):] = [np.zeros(len(x.columns))]

for j in range(N) :
    cv = StratifiedKFold(n_splits = K, shuffle = True)    
    for i, x in enumerate(tx) :
        x = x.to_numpy()
        for k, (train, test) in enumerate(cv.split(x, y)) :
            if(i % 2):
                xS, yS = sme.fit_sample(x[train], y[train])
            else:
                xS, yS = x[train], y[train]
            xtest, ytest = x[test], y[test]
            clf[i].fit(xS, yS)
            probas = clf[i].predict_proba(xtest)
            fpr, tpr, _ = roc_curve(ytest, probas[:, 1]) 
            mean_tprK[j, :, i] += np.interp(mean_fpr, fpr, tpr)
            mean_tprK[j, 0, i] = 0.0
            intercept[j, i] += clf[i].intercept_
            coef[i][j, :] += clf[i].coef_.ravel()
            
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

#%%
fig, ax_roc = plt.subplots() 
tcolor = [darkPurple, darkBlue, orangePurple, orangePurple, orange, orange,
          greenBlue, greenBlue]

for i in [0, 1, 2, 4, 6] :
    ax_roc.plot(mean_fpr, mean_tprN[:, i], linewidth = 6, color = tcolor[i])
    ci = stats.t.interval(level, dof, mean_tprN[:, i], std_tprN[:, i] / np.sqrt(N))
    ax_roc.fill_between(mean_fpr, ci[0], ci[1], color = tcolor[i],
                        alpha = 0.1)
    
ax_roc.plot([0, 1], [0, 1], '--', color = 'tab:gray', linewidth = 6)
ax_roc.set(xlabel = 'FPR', ylabel = 'TPR', ylim = [0, 1])
mean_auc_roc = np.mean(auc_roc, axis = 0)
std_auc_roc = np.std(auc_roc, axis = 0)

ax_roc.legend(['Prediction 1 (AUC = %0.2f $\pm$ %0.2f)' 
               % (mean_auc_roc[0], std_auc_roc[0]), 
               'Prediction 1* (AUC = %0.2f $\pm$ %0.2f)' 
               % (mean_auc_roc[1], std_auc_roc[1]),
               'Prediction 2 (AUC = %0.2f $\pm$ %0.2f)' 
               % (mean_auc_roc[2], std_auc_roc[2]),
               'Prediction 3 (AUC = %0.2f $\pm$ %0.2f)' 
               % (mean_auc_roc[4], std_auc_roc[4]), 
               'Prediction 4 (AUC = %0.2f $\pm$ %0.2f)' 
               % (mean_auc_roc[6], std_auc_roc[6])], loc = 'lower right',
               fontsize = 21)

#%%
fig, ax = plt.subplots() 

ax = sns.boxplot(data = scoresMean, orient = 'v',
                 linewidth = 2, fliersize = 6)
#ax = sns.swarmplot(data = scoresMean, color = ".25")
add_stat_annotation(ax, data = scoresMean, box_pairs = [('im_feat.','im_feat_OS'),
                                                  ('tissue_feat.', 'tissue_feat_OS'),
                                                  ('init_tum_area', 'init_tum_area_OS'),
                                                  ('8w_tum_area_norm', '8w_tum_area_OS_norm')],
                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
                    linewidth = 2, line_offset = -0.05, verbose = 2)

ax.set(ylabel = 'AUC')
ax.set_xticklabels([], ha = 'center')
#%%
plt.rcParams.update({'font.size': 32})
N = 100
K = 3
logReg = LogisticRegression()

#tx = [data[['tum_vol', 'ADC_ave', 'T2w_ave']],
#      data[['tum_area_from_vol', 'dens_ADCT2w']],
#      data[['init_tum_area']],
#      data[['TTum85']]]

x = data[['tum_vol', 'ADC_ave', 'T2w_ave']]
#x = data[['TTum330_alphaG1120']]
y = data[['bio_rec']].to_numpy().ravel() 


tpipeline = [Pipeline(steps = [('model', logReg)]),
             Pipeline(steps = [('over', SMOTE()), ('model', logReg)]),
             Pipeline(steps = [('over', BorderlineSMOTE()), ('model', logReg)]),
             Pipeline(steps = [('over', SMOTEENN()), ('model', logReg)]),
             Pipeline(steps = [('over', ADASYN()), ('model', logReg)])] 

scores = np.zeros((N, K, len(tpipeline)))
    
for j in range(N) :
        cv = StratifiedKFold(n_splits = K, shuffle = True)
        for i, pipeline in enumerate(tpipeline) :
            scores[j, :, i] = cross_val_score(pipeline, x, y, cv = cv,
                  scoring = 'roc_auc')
#            scores[j, :, i] = cross_val_score(logReg, x, y, cv = cv,
#                  scoring = 'roc_auc')
        
scoresMean = np.mean(scores, axis = 1)
#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['No oversample', 'SMOTE', 'BD SMOTE',
#                                     'SMOTE ENN', 'ADASYN'])
    
scoresMean = pd.DataFrame(data = scoresMean,
                          columns = ['No oversample', '04', '05', '06'])

print('Mean AUC')
print(scoresMean.mean(axis = 0))

print('Standard deviation AUC')
print(scoresMean.std(axis = 0))

print('Median AUC')
print(scoresMean.median(axis = 0))

#%%
plt.rcParams.update({'font.size': 32})
N = 100
K = 3
logReg = LogisticRegression()

#tx = [data[['tum_vol', 'ADC_ave', 'T2w_ave']],
#      data[['tum_area_from_vol', 'dens_ADCT2w']],
#      data[['init_tum_area']],
#      data[['TTum85']]]

x = data[['tum_vol', 'ADC_ave', 'T2w_ave']]
#x = data[['TTum330_alphaG1120']]
y = data[['bio_rec']].to_numpy().ravel() 


tpar = ['all'] 
#tpar = [3, 4, 5] 

scores = np.zeros((N, K, len(tpar) + 1))
    
for j in range(N) :
        cv = StratifiedKFold(n_splits = K, shuffle = True)
        scores[j, :, 0] = cross_val_score(logReg, x, y, cv = cv,
              scoring = 'roc_auc')
        for i, par in enumerate(tpar) :
            pipeline = Pipeline(steps = [('over', SMOTEENN(sampling_strategy = par)),
                                         ('model', logReg)])
            scores[j, :, i + 1] = cross_val_score(pipeline, x, y, cv = cv,
                  scoring = 'roc_auc')

        
scoresMean = np.mean(scores, axis = 1)   
scoresMean = pd.DataFrame(data = scoresMean,
                          columns = ['No oversample', 'all'])

print('Mean AUC')
print(scoresMean.mean(axis = 0))

print('Standard deviation AUC')
print(scoresMean.std(axis = 0))

print('Median AUC')
print(scoresMean.median(axis = 0))