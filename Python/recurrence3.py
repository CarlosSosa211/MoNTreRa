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
from sklearn.metrics import auc, confusion_matrix, roc_curve
from sklearn.metrics import precision_recall_curve, classification_report
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold, StratifiedKFold, LeaveOneOut
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.calibration import calibration_curve
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict 
from imblearn.over_sampling import SMOTE
from imblearn.combine import SMOTEENN

#%%
def fsigmoid(x, a, b):
    return 1.0 / (1.0 + np.exp(-a *(x - b)))
#    return np.exp(-a * np.exp(-b * x + c))
    
#%%
green = [153/255, 255/255, 102/255]
darkGreen = [127/255, 207/255, 127/255]
lightGreen = [205/255, 255/255, 105/255]
blue = [127/255, 185/255, 225/255]
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
#data = pd.read_csv('../../Carlos/Results/Recurrence/noHypNec_vascDensNoPref0.05/'
#                   'noHypNec_vascDensNoPref_05.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/vascDensUniPref0.03_ADC/'
#                   'vascDensUniPref_03.csv')
data = pd.read_csv('../../Carlos/Results/Recurrence/rec_summary_8wTumVol.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/simp/rec_summary_8wTumVol.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/simp/rec_summary_8wIntTumVol.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120/TTum330_alphaG1120.csv')

#%%
plt.close('all')
nfig = 0

fig, ax = plt.subplots()
sns.displot(data[['vascDensUniPref_03_ADC']].to_numpy())
sns.displot(data[['vascDensUniPref_03_simp_ADC']].to_numpy())
ax.set(title = 'Density', xlabel = '',
       ylabel = 'density')
ax.legend(['vascDensUniPref_03_ADC', 'vascDensUniPref_03_simp_ADC'])

_, p = stats.ttest_ind(data[['vascDensUniPref_03_ADC']].to_numpy(),
                      data[['vascDensUniPref_03_simp_ADC']].to_numpy())
print('vascDensUniPref_03', p)

fig, ax = plt.subplots()
sns.displot(data[['vascDensUniPref_038_ADC']].to_numpy())
sns.displot(data[['vascDensUniPref_038_simp_ADC']].to_numpy())
ax.set(title = 'Density', xlabel = '',
       ylabel = 'density')
ax.legend(['vascDensUniPref_038_ADC', 'vascDensUniPref_038_simp_ADC'])

_, p = stats.ttest_ind(data[['vascDensUniPref_038_ADC']].to_numpy(),
                      data[['vascDensUniPref_038_simp_ADC']].to_numpy())
print('vascDensUniPref_038', p)

fig, ax = plt.subplots()
sns.displot(data[['vascDensNoPref_05_ADC']].to_numpy())
sns.displot(data[['vascDensNoPref_05_simp_ADC']].to_numpy())
ax.set(title = 'Density', xlabel = '',
       ylabel = 'density')
ax.legend(['vascDensNoPref_05_ADC', 'vascDensNoPref_05_simp_ADC'])

_, p = stats.ttest_ind(data[['vascDensNoPref_05_ADC']].to_numpy(),
                      data[['vascDensNoPref_05_simp_ADC']].to_numpy())
print('vascDensNoPref_05', p)

fig, ax = plt.subplots()
sns.displot(data[['noHypNec_vascDensNoPref_03_ADC']].to_numpy())
sns.displot(data[['noHypNec_vascDensNoPref_03_simp_ADC']].to_numpy())
ax.set(title = 'Density', xlabel = '',
       ylabel = 'density')
ax.legend(['noHypNec_vascDensNoPref_03_ADC',
           'noHypNec_vascDensUniPref_03_simp_ADC'])

_, p = stats.ttest_ind(data[['noHypNec_vascDensNoPref_03_ADC']].to_numpy(),
                      data[['noHypNec_vascDensNoPref_03_simp_ADC']].to_numpy())
print('noHypNec_vascDensNoPref_03', p)

fig, ax = plt.subplots()
sns.displot(data[['noHypNec_vascDensNoPref_038_ADCT2w']].to_numpy())
sns.displot(data[['noHypNec_vascDensNoPref_038_simp3_ADCT2w']].to_numpy())
ax.set(title = 'Density', xlabel = '',
       ylabel = 'density')
ax.legend(['noHypNec_vascDensNoPref_038_ADCT2w',
           'noHypNec_vascDensUniPref_038_simp3_ADCT2w'])

_, p = stats.ttest_ind(data[['noHypNec_vascDensNoPref_038_ADCT2w']].to_numpy(),
                      data[['noHypNec_vascDensNoPref_038_simp3_ADCT2w']].to_numpy())
print('noHypNec_vascDensNoPref_038', p)


#%%
plt.close('all')
nfig = 0

ax = sns.displot(data = data, x = 'T2w_diff_var_mean', y = 'T2w_contrast_mean', hue = 'bio_rec', kind = 'kde')
ax.set_axis_labels('T2-w diff. var. mean', 'T2-w contrast mean')

ax = sns.displot(data = data, x = 'tum_vol', y = 'T2w_contrast_mean', hue = 'bio_rec', kind = 'kde')
ax.set_axis_labels('Tumor volume', 'T2-w contrast mean')

ax = sns.displot(data = data, x = 'ADC_ave', y = 'T2w_ave', hue = 'bio_rec', kind = 'kde')
ax.set_axis_labels('Average ADC', 'Average T2-w')

ax = sns.displot(data = data, x = 'tum_area_from_vol', y = 'dens_ADCT2w', hue = 'bio_rec', kind = 'kde')
ax.set_axis_labels('Tumor area (mm²)', 'Cell density from T2-w and ADC')

ax = sns.displot(data = data, x = 'init_tum_area', hue = 'bio_rec', kde = True)
ax.set_axis_labels('Initial tumor area', 'Density')
_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['init_tum_area'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['init_tum_area'].to_numpy())
print('init_tum_area', p)

ax = sns.displot(data = data, x = 'killed_90', hue = 'bio_rec', kde = True)
ax.set_axis_labels('Killed 90', 'Density')
_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['killed_90'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['killed_90'].to_numpy())
print('killed_90', p)

ax = sns.displot(data = data, x = 'TTum330_alphaG1120', hue = 'bio_rec', kde = True)
_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['TTum330_alphaG1120'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['TTum330_alphaG1120'].to_numpy())
print('TTum330_alphaG1120', p)

ax = sns.displot(data = data, x = 'TTum330_alphaG1120_norm', hue = 'bio_rec', kde = True)
_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['TTum330_alphaG1120_norm'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['TTum330_alphaG1120_norm'].to_numpy())
print('TTum330_alphaG1120_norm', p)

ax = sns.displot(data = data, x = 'TTum330_alphaG1120_immuno0.1_5', hue = 'bio_rec', kde = True)
_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['TTum330_alphaG1120_immuno0.1_5'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['TTum330_alphaG1120_immuno0.1_5'].to_numpy())
print('TTum330_alphaG1120_immuno0.1_5', p)

ax = sns.displot(data = data, x = 'TTum330_alphaG1120_immuno0.1_5_norm', hue = 'bio_rec', kde = True)
_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['TTum330_alphaG1120_immuno0.1_5_norm'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['TTum330_alphaG1120_immuno0.1_5_norm'].to_numpy())
print('TTum330_alphaG1120_immuno0.1_5_norm', p)



#%%
fig, ax = plt.subplots()
sns.displot(data.loc[data['bio_rec'] == 1]['init_tum_area'].to_numpy())
sns.displot(data.loc[data['bio_rec'] == 0]['init_tum_area'].to_numpy())
fig, ax = plt.subplots()
sns.displot(data.loc[data['bio_rec'] == 1]['8w_int_tum_area_norm'].to_numpy())
sns.displot(data.loc[data['bio_rec'] == 0]['8w_int_tum_area_norm'].to_numpy())

ax.legend(['Bio. rec. init. tum. area', 'No bio. rec. init. tum. area',
           'Bio. rec. 8w int. tum. area norm.', 'No bio. rec. 8w int. tum. area norm.'])



#%%
fig, ax = plt.subplots()
plt.scatter(data.loc[data['bio_rec'] == 0]['T2w_diff_var_mean'].to_numpy(),
            data.loc[data['bio_rec'] == 0]['T2w_contrast_mean'].to_numpy())
plt.scatter(data.loc[data['bio_rec'] == 1]['T2w_diff_var_mean'].to_numpy(),
            data.loc[data['bio_rec'] == 1]['T2w_contrast_mean'].to_numpy())

ax.set(xlabel = 'T2w_diff_var_mean', ylabel = 'T2w_contrast_mean')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

fig, ax = plt.subplots()
plt.scatter(data.loc[data['bio_rec'] == 0]['tum_vol'].to_numpy(),
            data.loc[data['bio_rec'] == 0]['T2w_contrast_mean'].to_numpy())
plt.scatter(data.loc[data['bio_rec'] == 1]['tum_vol'].to_numpy(),
            data.loc[data['bio_rec'] == 1]['T2w_contrast_mean'].to_numpy())

ax.set(xlabel = 'tum_vol', ylabel = 'T2w_contrast_mean')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

fig, ax = plt.subplots()
plt.scatter(data.loc[data['bio_rec'] == 0]['tum_area_from_vol'].to_numpy(),
            data.loc[data['bio_rec'] == 0]['dens_ADCT2w'].to_numpy())
plt.scatter(data.loc[data['bio_rec'] == 1]['tum_area_from_vol'].to_numpy(),
            data.loc[data['bio_rec'] == 1]['dens_ADCT2w'].to_numpy())

ax.set(title = 'Density', xlabel = 'tum_area_from_vol', ylabel = 'dens_ADCT2w')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

fig, ax = plt.subplots()
plt.scatter(data.loc[data['bio_rec'] == 0]['ADC_ave'].to_numpy(),
            data.loc[data['bio_rec'] == 0]['dens_ADC'].to_numpy())
plt.scatter(data.loc[data['bio_rec'] == 1]['ADC_ave'].to_numpy(),
            data.loc[data['bio_rec'] == 1]['dens_ADC'].to_numpy())

ax.set(title = 'Density', xlabel = 'ADC_ave', ylabel = 'T2w_ave')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

fig, ax = plt.subplots()
plt.scatter(data.loc[data['bio_rec'] == 0]['dens_ADC'].to_numpy(),
            data.loc[data['bio_rec'] == 0]['dens_T2w'].to_numpy())
plt.scatter(data.loc[data['bio_rec'] == 1]['dens_ADC'].to_numpy(),
            data.loc[data['bio_rec'] == 1]['dens_T2w'].to_numpy())

ax.set(title = 'Density', xlabel = 'dens_ADC', ylabel = 'dens_T2w')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

fig, ax = plt.subplots()
plt.scatter(data.loc[data['bio_rec'] == 0]['init_tum_area'].to_numpy(),
            data.loc[data['bio_rec'] == 0]['init_tum_area'].to_numpy())
plt.scatter(data.loc[data['bio_rec'] == 1]['init_tum_area'].to_numpy(),
            data.loc[data['bio_rec'] == 1]['init_tum_area'].to_numpy())

ax.set(title = 'Density', xlabel = 'init_tum_area',
       ylabel = 'init_tum_area', xlim = [0, 8], ylim = [0, 8])
ax.legend(['No. bio. rec.', 'Bio. rec.'])

fig, ax = plt.subplots()
plt.scatter(data.loc[data['bio_rec'] == 0]['8w_tum_area'].to_numpy(),
            data.loc[data['bio_rec'] == 0]['8w_tum_area'].to_numpy())
plt.scatter(data.loc[data['bio_rec'] == 1]['8w_tum_area'].to_numpy(),
            data.loc[data['bio_rec'] == 1]['8w_tum_area'].to_numpy())

ax.set(title = 'Density', xlabel = '8w_tum_area',
       ylabel = '8w_tum_area', xlim = [0, 8], ylim = [0, 8])
ax.legend(['No. bio. rec.', 'Bio. rec.'])

fig, ax = plt.subplots()
plt.scatter(data.loc[data['bio_rec'] == 0]['TTum330_alphaG1120_immuno0.1_5_norm'].to_numpy(),
            data.loc[data['bio_rec'] == 0]['TTum330_alphaG1120_immuno0.1_5_norm'].to_numpy())
plt.scatter(data.loc[data['bio_rec'] == 1]['TTum330_alphaG1120_immuno0.1_5_norm'].to_numpy(),
            data.loc[data['bio_rec'] == 1]['TTum330_alphaG1120_immuno0.1_5_norm'].to_numpy())

ax.set(title = 'Density', xlabel = 'TTum330_alphaG1120_immuno0.1_5_norm',
       ylabel = 'TTum330_alphaG1120_immuno0.1_5_norm')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

fig, ax = plt.subplots()
plt.scatter(data.loc[data['bio_rec'] == 0]['8w_int_tum_area_norm'].to_numpy(),
            data.loc[data['bio_rec'] == 0]['8w_int_tum_area_norm'].to_numpy())
plt.scatter(data.loc[data['bio_rec'] == 1]['8w_int_tum_area_norm'].to_numpy(),
            data.loc[data['bio_rec'] == 1]['8w_int_tum_area_norm'].to_numpy())

ax.set(title = 'Density', xlabel = '8w_int_tum_area_norm',
       ylabel = '8w_int_tum_area_norm')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

fig, ax = plt.subplots()
plt.scatter(data.loc[data['bio_rec'] == 0]['killed_90'].to_numpy(),
            data.loc[data['bio_rec'] == 0]['killed_90'].to_numpy())
plt.scatter(data.loc[data['bio_rec'] == 1]['killed_90'].to_numpy(),
            data.loc[data['bio_rec'] == 1]['killed_90'].to_numpy())

ax.set(title = 'Density', xlabel = 'killed_90',
       ylabel = 'killed_90')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

#%%
plt.close('all')
nfig = 0

fig, ax = plt.subplots()
sns.displot(data.loc[data['bio_rec'] == 1]['ADC_ave'].to_numpy())
sns.displot(data.loc[data['bio_rec'] == 0]['ADC_ave'].to_numpy())
ax.set(title = 'Density', xlabel = 'Initial ADC', ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['ADC_ave'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['ADC_ave'].to_numpy())
print('Initial ADC', p)

fig, ax = plt.subplots()
sns.displot(data.loc[data['bio_rec'] == 1]['max_tum_area'].to_numpy())
sns.displot(data.loc[data['bio_rec'] == 0]['max_tum_area'].to_numpy())
ax.set(title = 'Density', xlabel = 'Initial maximal tumor area',
       ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['max_tum_area'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['max_tum_area'].to_numpy())
print('Initial maximal tumor area', p)

fig, ax = plt.subplots()
sns.displot(data.loc[data['bio_rec'] == 1]['psa'].to_numpy())
sns.displot(data.loc[data['bio_rec'] == 0]['psa'].to_numpy())
ax.set(title = 'Density', xlabel = 'Initial maximal tumor area',
       ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['psa'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['psa'].to_numpy())
print('Initial PSA', p)

fig, ax = plt.subplots()
sns.displot(data.loc[data['bio_rec'] == 1]['two_mon_tum_area'].to_numpy())
sns.displot(data.loc[data['bio_rec'] == 0]['two_mon_tum_area'].to_numpy())
ax.set(title = 'Density', xlabel = '$In$ $silico$ 2 months tumor area',
       ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['two_mon_tum_area'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['two_mon_tum_area'].to_numpy())
print('In silico 2 months tumor area', p)
    
#%%

data.boxplot(by = 'bio_rec_6', column = ['8-500'], grid = False,
             fontsize = 32)

data.boxplot(by = 'bio_rec_6', column = ['8-500'], grid = False,
             fontsize = 32)
#%%
nfig += 1
plt.figure(nfig)
data.boxplot(by = 'bio_rec', column = ['max_tum_area'], grid = False,
             fontsize = 32)


#%%
plt.close('all')
fig, ax = plt.subplots()
adcAve = data[['ADC_ave']].to_numpy().ravel() 
T2wAve = data[['T2w_ave']].to_numpy().ravel()
ax.plot(adcAve, T2wAve, 'o', linewidth = 6)
p = np.polyfit(adcAve, T2wAve, 1)
ax.plot(adcAve, np.polyval(p, adcAve), '--k', linewidth = 6)
ax.set(xlabel = 'average ADC', ylabel = 'average T2w')
    
fig, ax = plt.subplots()    
ax.plot(data[['epi_nuc_dens_norm_ADC_ave']],
        data[['epi_nuc_dens_norm_T2w_ave']], 'o',
        linewidth = 6)
ax.plot([-1, 1], [-1, 1], '--k', linewidth = 6)
ax.set(xlabel = 'normalized epithelial nucleus density (from  average ADC)',
       ylabel = 'normalized epithelial nucleus density (from average T2w)')
    
#%%
fig, ax = plt.subplots()    
ax.plot(data[['ADC_average']],
        data[['tumor_volume']], 'o',
        linewidth = 6)
ax.set(xlabel = 'average ADC',
       ylabel = 'tumor volume')

#%%
pca = PCA(n_components = 10)
pca.fit(X)
X = pca.fit_transform(X)

nfig += 1
plt.figure(nfig)
plt.plot(pca.explained_variance_ratio_)
plt.ylabel("Explained variance")
plt.xlabel("Number of components")
plt.show()

#%%
plt.close('all')

data = np.loadtxt('../InputFiles/tumDensADCT2w.dat')
nrow, _ = np.shape(data)

yfromADC = data[:, 2]
yfromT2w = data[:, 3]
ymean = 0.5 * (yfromADC + yfromT2w)


xfromADC = np.zeros([nrow, 2])
xfromADC[:, 0] = data[:, 0]
xfromADC[:, 1] = yfromADC / -0.39

xfromT2w = np.zeros([nrow, 2])
xfromT2w[:, 0] = yfromT2w / -0.56
xfromT2w[:, 1] = data[:, 1]

xmean = np.zeros([nrow, 2])
xmean[:, 0] = ymean / -0.56
xmean[:, 1] = ymean / -0.39


ax = plt.axes(projection = '3d')
ax.plot3D(xfromADC[:, 0].ravel(), xfromADC[:,1].ravel(), yfromADC.ravel(), 'o',
          linewidth = 4)
ax.plot3D(xfromT2w[:, 0].ravel(), xfromT2w[:,1].ravel(), yfromT2w.ravel(), 'o',
          linewidth = 4)
ax.plot3D(xmean[:, 0].ravel(), xmean[:,1].ravel(), ymean.ravel(), 'o',
          linewidth = 4)
ax.set(title = 'Epithelial densitiy', xlabel = 'ADC', ylabel = 'T2w')
#ax.legend(['Data', 'Prediction'])

#x = np.concatenate((xfromADC, xfromT2w))
#y = np.concatenate((yfromADC, yfromT2w))
#interp = LinearNDInterpolator(x, y)
#
#yinterp = interp._call_(x)

