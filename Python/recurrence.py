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
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import StratifiedKFold
from imblearn.over_sampling import SMOTE

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
data = pd.read_csv('../../Carlos/Results/Recurrence/simp/rec_summary_8wTumVol.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/simp/rec_summary_8wIntTumVol.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/rec_summary_8wTumDens.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/rec_summary_8wIntTumDens.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/rec_summary_12wIntTumVol.csv')

#%%
plt.rcParams.update({'font.size': 32})
N = 1000
K = 3
#logReg = LogisticRegression()
logReg = MLPClassifier(activation = 'logistic', max_iter = 1000)
#data[['ADC_med', 'max_tum_area', 'T2w_diff_var_mean', 'tum_vol',
#            'T2w_contrast_mean']]

tx = [data[['tum_vol', 'ADC_ave', 'T2w_ave', 'total_dose']],
      data[['tum_area_from_vol', 'dens_ADCT2w', 'total_dose']],
      data[['init_tum_area_tumVolADCT2w', 'total_dose']],
      data[['TTum330_alphaG1120_ADCT2w_norm']]]

#tx = [data[['TTum260_alphaG1120_ADCT2w']],
#      data[['TTum330_alphaG1120_ADCT2w']],
#      data[['TTum360_alphaG1120_ADCT2w']],
#      data[['TTum400_alphaG1120_ADCT2w']]]

#tx = [data[['TTum330_alphaG1090_ADCT2w']],
#      data[['TTum330_alphaG1120_ADCT2w']],
#      data[['TTum330_ADCT2w']],
#      data[['TTum330_alphaG1223_ADCT2w']]]

y = data[['bio_rec']].to_numpy().ravel()  

scores = np.zeros((N, K, len(tx)))
    
for j in range(N) :
        cv = StratifiedKFold(n_splits = K, shuffle = True)
        for i, x in enumerate(tx) :
            scores[j, :, i] = cross_val_score(logReg, x, y, cv = cv,
                  scoring = 'roc_auc')
        
scoresMean = np.mean(scores, axis = 1)
scoresMean = pd.DataFrame(data = scoresMean,
                          columns = ['im_feat.', 'tissue_feat.',
                                     'init_tum_area', '8w_tum_area_norm'])
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
plt.close('all')

tcolor = [red, orange, blue, green]
ax = sns.boxplot(data = scoresMean, orient = 'v', palette = tcolor)
#ax = sns.swarmplot(data = scoresMean, color = ".25")
add_stat_annotation(ax, data = scoresMean,
                    box_pairs = [('im_feat.',
                                'tissue_feat.'),
                                ('tissue_feat.',
                                'init_tum_area'),
                                ('init_tum_area',
                                '8w_tum_area_norm')],
                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
                    line_offset = -0.05, verbose = 2)
#add_stat_annotation(ax, data = scoresMean,
#                    box_pairs = [('260',
#                                '330'),
#                                ('330',
#                                '400'),
#                                ('260',
#                                '400')],
#                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
#                    line_offset = -0.05, verbose = 2)
ax.set(ylabel = 'AUC', ylim = [0.0, 1])
ax.set_xticklabels([], ha = 'center')

#%%   
plt.close('all')
plt.rcParams.update({'font.size': 32})  
N = 100
K = 3

#tx = [data[['tum_vol', 'ADC_ave', 'T2w_ave']],
#      data[['tum_area_from_vol', 'dens_ADCT2w']],
#      data[['init_tum_area_tumVolADCT2w']],
#      data[['TTum330_alphaG1120_ADCT2w']]]

tx = [data[['tum_vol', 'ADC_ave', 'T2w_ave']],
      data[['tum_area_from_vol', 'dens_ADCT2w']],
      data[['init_tum_area_tumVolADCT2w']],
      data[['TTum330_alphaG1120_ADCT2w_norm']]]

y = data[['bio_rec']].to_numpy().ravel()  

mean_fpr = np.linspace(0, 1, 500)
mean_tprK = np.zeros((N, 500, len(tx)))
auc_ = np.zeros((N, len(tx))) 

logReg = LogisticRegression()

for j in range(N) :
    cv = StratifiedKFold(n_splits = K, shuffle = True)    
    for i, x in enumerate(tx) :
        x = x.to_numpy()
        for k, (train, test) in enumerate(cv.split(x, y)) :
            xS, yS = SMOTE().fit_sample(x[train], y[train])
#            xS, yS = x[train], y[train]
            logReg.fit(xS, yS)
            probas = logReg.predict_proba(x[test])
            fpr, tpr, _ = roc_curve(y[test], probas[:, 1]) 
            mean_tprK[j, :, i] += np.interp(mean_fpr, fpr, tpr)
            mean_tprK[j, 0, i] = 0.0
            
mean_tprK /= K
mean_tprK[:, -1, :] = 1.0
auc_ = np.trapz(mean_tprK, mean_fpr, axis = 1) 
    
mean_tprN = np.mean(mean_tprK, axis = 0)
std_tprN = np.std(mean_tprK, axis = 0)
level = 0.95
dof = K - 1

fig, ax = plt.subplots() 
tcolor = [red, orange, blue, green]

for i in range(len(tx)) :
    ax.plot(mean_fpr, mean_tprN[:, i], linewidth = 6, color = tcolor[i])
    ci = stats.t.interval(level, dof, mean_tprN[:, i], std_tprN[:, i])
    ax.fill_between(mean_fpr, ci[0], ci[1], color = tcolor[i],
                    alpha = 0.1)

ax.plot([0, 1], [0, 1], '--', color = 'tab:gray', linewidth = 6)
ax.set(xlabel = 'FPR', ylabel = 'TPR', ylim = [0, 1])
mean_auc = np.mean(auc_, axis = 0)
std_auc = np.std(auc_, axis = 0)
#ax.legend(["Pre-treatment\nimaging parameters\n"
#           "(AUC = %0.2f $\pm$ %0.2f)" % (mean_auc[0], std_auc[0]),
#           "Tum. area at 8 w. from\ncomprehensive model\n"
#           "(AUC = %0.2f $\pm$ %0.2f)" % (mean_auc[1], std_auc[1]),
#           "Tum. area at 8 w. from\nreduced model\n"
#           "(AUC = %0.2f $\pm$ %0.2f)" % (mean_auc[2], std_auc[2])],
##           '(AUC = %0.2f)' % mean_auc[3],
##           '(AUC = %0.2f)' % mean_auc[4],
##           '(AUC = %0.2f)' % mean_auc[5],
##           '(AUC = %0.2f)' % mean_auc[6]],
#           loc = 'lower right')
#
#           "Tum. area at 0 w.\n"
#           "(AUC = %0.2f $\pm$ %0.2f)" % (mean_auc[1], std_auc[1]),
ax.legend(["Imaging parameters\n(AUC = %0.2f $\pm$ %0.2f)" %
           (mean_auc[0], std_auc[0]),
           "Tissue parameters\n(AUC = %0.2f $\pm$ %0.2f)" %
           (mean_auc[1], std_auc[1]),
           "Tum. area at 0 w.\n(AUC = %0.2f $\pm$ %0.2f)" % 
           (mean_auc[2], std_auc[2]),
           "Tum. area at 8 w.\n(AUC = %0.2f $\pm$ %0.2f)" %
           (mean_auc[3], std_auc[3])],
           loc = 'lower right')
#
#_, p = stats.wilcoxon(auc_[:, 0], auc_[:, 1], alternative = 'less')
#print("4_features, comp.:", p)
#_, p = stats.wilcoxon(auc_[:, 0], auc_[:, 2], alternative = 'less')
#print("4_features, red.:", p)
#_, p = stats.wilcoxon(auc_[:, 1], auc_[:, 2])
#print("comp., red.:", p)

#%%
plt.close('all')
nfig = 0

fig, ax = plt.subplots()
sns.distplot(data[['vascDensUniPref_03_ADC']].to_numpy())
sns.distplot(data[['vascDensUniPref_03_simp_ADC']].to_numpy())
ax.set(title = 'Density', xlabel = '',
       ylabel = 'density')
ax.legend(['vascDensUniPref_03_ADC', 'vascDensUniPref_03_simp_ADC'])

_, p = stats.ttest_ind(data[['vascDensUniPref_03_ADC']].to_numpy(),
                      data[['vascDensUniPref_03_simp_ADC']].to_numpy())
print('vascDensUniPref_03', p)

fig, ax = plt.subplots()
sns.distplot(data[['vascDensUniPref_038_ADC']].to_numpy())
sns.distplot(data[['vascDensUniPref_038_simp_ADC']].to_numpy())
ax.set(title = 'Density', xlabel = '',
       ylabel = 'density')
ax.legend(['vascDensUniPref_038_ADC', 'vascDensUniPref_038_simp_ADC'])

_, p = stats.ttest_ind(data[['vascDensUniPref_038_ADC']].to_numpy(),
                      data[['vascDensUniPref_038_simp_ADC']].to_numpy())
print('vascDensUniPref_038', p)

fig, ax = plt.subplots()
sns.distplot(data[['vascDensNoPref_05_ADC']].to_numpy())
sns.distplot(data[['vascDensNoPref_05_simp_ADC']].to_numpy())
ax.set(title = 'Density', xlabel = '',
       ylabel = 'density')
ax.legend(['vascDensNoPref_05_ADC', 'vascDensNoPref_05_simp_ADC'])

_, p = stats.ttest_ind(data[['vascDensNoPref_05_ADC']].to_numpy(),
                      data[['vascDensNoPref_05_simp_ADC']].to_numpy())
print('vascDensNoPref_05', p)

fig, ax = plt.subplots()
sns.distplot(data[['noHypNec_vascDensNoPref_03_ADC']].to_numpy())
sns.distplot(data[['noHypNec_vascDensNoPref_03_simp_ADC']].to_numpy())
ax.set(title = 'Density', xlabel = '',
       ylabel = 'density')
ax.legend(['noHypNec_vascDensNoPref_03_ADC',
           'noHypNec_vascDensUniPref_03_simp_ADC'])

_, p = stats.ttest_ind(data[['noHypNec_vascDensNoPref_03_ADC']].to_numpy(),
                      data[['noHypNec_vascDensNoPref_03_simp_ADC']].to_numpy())
print('noHypNec_vascDensNoPref_03', p)

fig, ax = plt.subplots()
sns.distplot(data[['noHypNec_vascDensNoPref_038_ADCT2w']].to_numpy())
sns.distplot(data[['noHypNec_vascDensNoPref_038_simp_ADCT2w']].to_numpy())
ax.set(title = 'Density', xlabel = '',
       ylabel = 'density')
ax.legend(['noHypNec_vascDensNoPref_038_ADCT2w',
           'noHypNec_vascDensUniPref_038_simp_ADCT2w'])

_, p = stats.ttest_ind(data[['noHypNec_vascDensNoPref_038_ADCT2w']].to_numpy(),
                      data[['noHypNec_vascDensNoPref_038_simp_ADCT2w']].to_numpy())
print('noHypNec_vascDensNoPref_05', p)


#%%
plt.close('all')
nfig = 0

fig, ax = plt.subplots()
sns.distplot(data.loc[data['bio_rec'] == 1]['ADC_ave'].to_numpy())
sns.distplot(data.loc[data['bio_rec'] == 0]['ADC_ave'].to_numpy())
ax.set(title = 'Density', xlabel = 'Initial ADC', ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['ADC_ave'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['ADC_ave'].to_numpy())
print('Initial ADC', p)

fig, ax = plt.subplots()
sns.distplot(data.loc[data['bio_rec'] == 1]['max_tum_area'].to_numpy())
sns.distplot(data.loc[data['bio_rec'] == 0]['max_tum_area'].to_numpy())
ax.set(title = 'Density', xlabel = 'Initial maximal tumor area',
       ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['max_tum_area'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['max_tum_area'].to_numpy())
print('Initial maximal tumor area', p)

fig, ax = plt.subplots()
sns.distplot(data.loc[data['bio_rec'] == 1]['vascDensUniPref_03'].to_numpy())
sns.distplot(data.loc[data['bio_rec'] == 0]['vascDensUniPref_03'].to_numpy())
ax.set(title = 'Density', xlabel = 'vascDensUniPref_03',
       ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['vascDensUniPref_03'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['vascDensUniPref_03'].to_numpy())
print('vascDensUniPref_03', p)

fig, ax = plt.subplots()
sns.distplot(data.loc[data['bio_rec'] == 1]['noHypNec_vascDensNoPref_03'].to_numpy())
sns.distplot(data.loc[data['bio_rec'] == 0]['noHypNec_vascDensNoPref_03'].to_numpy())
ax.set(title = 'Density', xlabel = 'noHypNec_vascDensNoPref_03',
       ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['noHypNec_vascDensNoPref_03'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['noHypNec_vascDensNoPref_03'].to_numpy())
print('noHypNec_vascDensNoPref_03', p)

fig, ax = plt.subplots()
sns.distplot(data.loc[data['bio_rec'] == 1]['noHypNec_vascDensNoPref_05'].to_numpy())
sns.distplot(data.loc[data['bio_rec'] == 0]['noHypNec_vascDensNoPref_05'].to_numpy())
ax.set(title = 'Density', xlabel = 'noHypNec_vascDensNoPref_05',
       ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['noHypNec_vascDensNoPref_05'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['noHypNec_vascDensNoPref_05'].to_numpy())
print('noHypNec_vascDensNoPref_05', p)

fig, ax = plt.subplots()
sns.distplot(data.loc[data['bio_rec'] == 1]['TTum330_alphaG1120_ADCT2w'].to_numpy())
sns.distplot(data.loc[data['bio_rec'] == 0]['TTum330_alphaG1120_ADCT2w'].to_numpy())
ax.set(title = 'Density', xlabel = '8w_tum_area',
       ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['vD038_tumVolADCT2w'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['vD038_tumVolADCT2w'].to_numpy())
print('vD038_tumVolADCT2w', p)

#%%
#fig, ax = plt.subplots()
#plt.scatter(data.loc[data['bio_rec'] == 0]['T2w_ave'].to_numpy(),
#         data.loc[data['bio_rec'] == 0]['vascDensUniPref_038_ADCT2w'].to_numpy())
#plt.scatter(data.loc[data['bio_rec'] == 1]['T2w_ave'].to_numpy(),
#         data.loc[data['bio_rec'] == 1]['vascDensUniPref_038_ADCT2w'].to_numpy())
#
#ax.set(title = 'Density', xlabel = 'ADC_ave', ylabel = 'T2w_ave')
#ax.legend(['No. bio. rec.', 'Bio. rec.'])

fig, ax = plt.subplots()
plt.scatter(data.loc[data['bio_rec'] == 0]['ADC_ave'].to_numpy(),
            data.loc[data['bio_rec'] == 0]['init_tum_area_ADC'].to_numpy())
plt.scatter(data.loc[data['bio_rec'] == 1]['ADC_ave'].to_numpy(),
            data.loc[data['bio_rec'] == 1]['init_tum_area_ADC'].to_numpy())

ax.set(title = 'Density', xlabel = 'ADC_ave', ylabel = 'init_tum_area_ADC')
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
plt.scatter(data.loc[data['bio_rec'] == 0]['TTum330_alphaG1120_ADCT2w'].to_numpy(),
            data.loc[data['bio_rec'] == 0]['TTum330_alphaG1120_ADCT2w'].to_numpy())
plt.scatter(data.loc[data['bio_rec'] == 1]['TTum330_alphaG1120_ADCT2w'].to_numpy(),
            data.loc[data['bio_rec'] == 1]['TTum330_alphaG1120_ADCT2w'].to_numpy())

ax.set(title = 'Density', xlabel = 'dens_ADC', ylabel = 'dens_T2w')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

#%%
plt.close('all')
nfig = 0

fig, ax = plt.subplots()
sns.distplot(data.loc[data['bio_rec'] == 1]['ADC_ave'].to_numpy())
sns.distplot(data.loc[data['bio_rec'] == 0]['ADC_ave'].to_numpy())
ax.set(title = 'Density', xlabel = 'Initial ADC', ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['ADC_ave'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['ADC_ave'].to_numpy())
print('Initial ADC', p)

fig, ax = plt.subplots()
sns.distplot(data.loc[data['bio_rec'] == 1]['max_tum_area'].to_numpy())
sns.distplot(data.loc[data['bio_rec'] == 0]['max_tum_area'].to_numpy())
ax.set(title = 'Density', xlabel = 'Initial maximal tumor area',
       ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['max_tum_area'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['max_tum_area'].to_numpy())
print('Initial maximal tumor area', p)

fig, ax = plt.subplots()
sns.distplot(data.loc[data['bio_rec'] == 1]['psa'].to_numpy())
sns.distplot(data.loc[data['bio_rec'] == 0]['psa'].to_numpy())
ax.set(title = 'Density', xlabel = 'Initial maximal tumor area',
       ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['psa'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['psa'].to_numpy())
print('Initial PSA', p)

fig, ax = plt.subplots()
sns.distplot(data.loc[data['bio_rec'] == 1]['two_mon_tum_area'].to_numpy())
sns.distplot(data.loc[data['bio_rec'] == 0]['two_mon_tum_area'].to_numpy())
ax.set(title = 'Density', xlabel = '$In$ $silico$ 2 months tumor area',
       ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['two_mon_tum_area'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['two_mon_tum_area'].to_numpy())
print('In silico 2 months tumor area', p)
    
#%%
plt.rcParams.update({'font.size': 32})
#tx = ['ADC_average', 'max_tumor_area', 'tumor_volume', 'T2w_average',
#      'epi_dens_ADC_ave', 'epi_dens_T2w_ave', 'epi_dens_mean']
#txNames = ['ADC', 'maxTumArea', 'tumVol', 'T2w', 'dens(ADC)', 'dens(T2w)',
#           'dens(ADC&T2w)']

#tx = ['ADC_average', 'max_tumor_area', 'T2w_average']
#txNames = ['average ADC', 'max_tumor_area', 'average T2w']
#
#tx = ['ADC_ave', 'ADC_ave_norm', 'epi_dens_norm_ADC_ave',
#      'epi_dens_ADC_ave', 'T2w_ave', 'T2w_ave_norm',
#      'epi_dens_norm_T2w_ave', 'epi_dens_T2w_ave', 'epi_dens_mean']
#txNames = ['ADC', 'ADCNorm', 'densNorm(ADC)', 'dens(ADC)', 'T2w', 'T2wNorm', 
#           'densNorm(T2w)', 'dens(T2w)', 'dens(ADC&T2w)']

#tx = ['ADC_ave', 'ADC_ave_norm', 'epi_dens_norm_ADC_ave',
#      'epi_dens_ADC_ave', 'T2w_ave', 'T2w_ave_norm',
#      'epi_dens_norm_T2w_ave', 'epi_dens_T2w_ave', 'epi_dens_mean']
#txNames = ['ADC', 'ADCNorm', 'densNorm(ADC)', 'dens(ADC)', 'T2w', 'T2wNorm', 
#           'densNorm(T2w)', 'dens(T2w)', 'dens(ADC&T2w)']
N = 1000
K = 3
logReg = LogisticRegression()

tx = [data[['max_tum_area', 'total_dose', 'ADC_ave']],
      data[['8w_tum_area']],
      data[['12w_tum_area']],
      data[['8w_int_tum_area']],
      data[['12w_int_tum_area']],
      data[['8w_tum_dens']],
      data[['12w_tum_dens']],
      data[['8w_int_tum_dens']],
      data[['12w_int_tum_dens']],
      data[['time_to_95']],
      data[['time_to_99']]]

#tx = [data[['max_tum_area', 'total_dose', 'ADC_ave']],
#      data[['max_tum_area', 'total_dose', 'ADC_ave', '8w_tum_area', '12w_tum_area',
#            '8w_int_tum_area', '12w_int_tum_area']]]
y = data[['bio_rec']].to_numpy().ravel() 

scores = np.zeros((N, K, len(tx)))
    
for j in range(N) :
        cv = StratifiedKFold(n_splits = K, shuffle = True)
        for i, x in enumerate(tx) :
            scores[j, :, i] = cross_val_score(logReg, x, y, cv = cv,
                  scoring = 'roc_auc')
        
scoresMean = np.mean(scores, axis = 1)
scoresMean = pd.DataFrame(data = scoresMean,
                          columns = ['3_features',
                                     '8w_tum_area',
                                     '12w_tum_area',
                                     '8w_int_tum_area',
                                     '12w_int_tum_area',
                                     '8w_tum_dens',
                                     '12w_tum_dens',
                                     '8w_int_tum_dens',
                                     '12w_int_tum_dens',
                                     'time_to_95',
                                     'time_to_99'])
    
scoresMean = pd.DataFrame(data = scoresMean, columns = ['3_features',
                                                        'in silico'])
    
#%%
plt.close('all')
ax = sns.boxplot(data = scoresMean, orient = 'v', palette = 'Set2')
ax = sns.swarmplot(data = scoresMean, color = ".25")
add_stat_annotation(ax, data = scoresMean,
                    box_pairs = [('3_features', 
                                 '8w_tum_area'),
                                ('3_features', 
                                 '12w_tum_area'),
                                ('3_features', 
                                 '8w_int_tum_area'),
                                ('3_features',
                                 '12w_int_tum_area'),
                                ('3_features', 
                                 '8w_tum_dens'),
                                ('3_features', 
                                 '12w_tum_dens'),
                                ('3_features', 
                                 '8w_int_tum_dens'),
                                ('3_features',
                                 '12w_int_tum_dens'),
                                ('3_features',
                                 'time_to_95'),
                                ('3_features',
                                 'time_to_99')],
                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
                    line_offset = -0.05, verbose = 2)
ax.set(title = 'AUC', ylabel = 'AUC', ylim = [None, 1])

#%%
nfig += 1
plt.figure(nfig)
data.boxplot(by = 'gleason', column = ['ADC_ave'], grid = False,
             fontsize = 32)

nfig += 1
plt.figure(nfig)
data.boxplot(by = 'gleason', column = ['tum_vol'], grid = False,
             fontsize = 32)

#%%
nfig += 1
plt.figure(nfig)
data.boxplot(by = 'bio_rec', column = ['max_tum_area'], grid = False,
             fontsize = 32)


#%%
plt.close('all')

sigmoid = np.zeros(301)
x = list(range(301))
fig, ax = plt.subplots();

adcRec = data.loc[data['bio_rec'] == 1][['TTum330_alphaG1120_ADCT2w']].to_numpy()
adcRec = adcRec.ravel()
ecdfRec = ECDF(adcRec)
popt, pcov = curve_fit(fsigmoid, ecdfRec.x[1:], ecdfRec.y[1:], p0 = [1, 50])
for j in range(301):
    sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
ax.plot(ecdfRec.x, ecdfRec.y, linewidth = 6)
#ax.plot(x, sigmoid, linewidth = 6)

adcNoRec = data.loc[data['bio_rec'] == 0][['TTum330_alphaG1120_ADCT2w']].to_numpy()
adcNoRec = adcNoRec.ravel()
ecdfNoRec = ECDF(adcNoRec)
popt, pcov = curve_fit(fsigmoid, ecdfNoRec.x[1:], ecdfNoRec.y[1:],
                       p0 = [1, 50])
for j in range(301):
    sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
ax.plot(ecdfNoRec.x, ecdfNoRec.y, linewidth = 6)
#ax.plot(x, sigmoid, linewidth = 6)
ax.set(xlabel = '8w_tum_area', ylabel = 'cumulative density', title = '8w_tum_area')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

thresholds = np.linspace(float(data[['TTum330_alphaG1120_ADCT2w']].min()),
                         float(data[['TTum330_alphaG1120_ADCT2w']].max()), 20)

fpr = np.zeros(len(thresholds))
tpr = np.zeros(len(thresholds))
for i, el in enumerate(thresholds) :
    fpr[i] = ecdfNoRec.y[ecdfNoRec.x < el][0]
    tpr[i] = ecdfRec.y[ecdfRec.x < el][0]

fig, ax = plt.subplots()
ax.plot(fpr, tpr, linewidth = 6)
ax.plot([0, 1], [0, 1], '--k', linewidth = 6)
ax.set(title = 'ROC')

#%%        
#    ax.plot(thresholds, tpr, linewidth = 6)
#    ax.set(title = 'Rates - ' + txNames[i], xlabel = 'threshold',
#           ylabel = 'rate')
#    ax.legend(['False positives', 'True positives'])

        
#    ypred = cross_val_predict(logReg, x, y, cv = 3)
    
    #xtrain, xtest, ytrain, ytest = train_test_split(x, y, stratify = y)
    
#    logReg = LogisticRegression()
#    logReg.fit(xtrain, ytrain)
#    ypred = logReg.predict(xtest)
#
#    confusion_matrix(ytest, ypred)
#
#    probas = logReg.predict_proba(xtest)
#    fpr, tpr, thresholds = roc_curve(ytest, probas[:, 0],
#                                 pos_label = logReg.classes_[0],
#                                 drop_intermediate = False)
#    rocAuc = auc(fpr, tpr)
#
#    fig, ax = plt.subplots()
#    ax.plot(thresholds, fpr, linewidth = 6)
#    ax.plot(thresholds, tpr, linewidth = 6)
#    ax.set(title = 'Rates - ' + txNames[i], xlabel = 'threshold',
#           ylabel = 'rate')
#    ax.legend(['False positives', 'True positives'])
#
#    fig, ax = plt.subplots()    
#    ax.plot(fpr, tpr, linewidth = 6)
#    ax.plot([0, 1], [0, 1], '--k', linewidth = 6)
#    ax.set(title = 'ROC - ' + txNames[i])
#    ax.legend(['ROC (AUC = %0.2f)' % rocAuc])



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

