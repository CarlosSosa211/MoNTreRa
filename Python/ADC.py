# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
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
from sklearn.linear_model import RidgeClassifier
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import StratifiedKFold
from imblearn.over_sampling import SMOTE
#from sklearn.feature_selection import VarianceThreshold
#from sklearn.decomposition import PCA
#from sklearn.model_selection import GridSearchCV
#from sklearn.linear_model import LinearRegression, Ridge, BayesianRidge, Lasso
#from sklearn.neighbors import KNeighborsClassifier
#from sklearn.tree import DecisionTreeRegressor
#from sklearn.ensemble import AdaBoostRegressor
#from sklearn.ensemble import RandomForestRegressor
#from sklearn.svm import SVR
#from sklearn.kernel_ridge import KernelRidge
#from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
#from sklearn.neural_network import MLPRegressor

#%%
#data = pd.read_csv('../../Carlos/Results/Recurrence/noHypNec_vascDensNoPref0.05/'
#                   'noHypNec_vascDensNoPref_05.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/vascDensUniPref0.03_ADC/'
#                   'vascDensUniPref_03.csv')
data = pd.read_csv('../../Carlos/Results/Recurrence/rec_summary_8wTumVol.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/rec_summary_8wIntTumVol.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/rec_summary_8wTumDens.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/rec_summary_8wIntTumDens.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/rec_summary_12wIntTumVol.csv')

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
sns.distplot(data.loc[data['bio_rec'] == 1]['noHypNec_vascDensNoPref_038_simp_ADCT2w'].to_numpy())
sns.distplot(data.loc[data['bio_rec'] == 0]['noHypNec_vascDensNoPref_038_simp_ADCT2w'].to_numpy())
ax.set(title = 'Density', xlabel = 'noHypNec_vascDensNoPref_038_simp',
       ylabel = 'density')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

_, p = stats.ttest_ind(data.loc[data['bio_rec'] == 1]['noHypNec_vascDensNoPref_038_simp_ADCT2w'].to_numpy(),
                      data.loc[data['bio_rec'] == 0]['noHypNec_vascDensNoPref_038_simp_ADCT2w'].to_numpy())
print('noHypNec_vascDensNoPref_038_simp', p)

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
plt.scatter(data.loc[data['bio_rec'] == 0]['noHypNec_vascDensNoPref_038_ADCT2w'].to_numpy(),
            np.zeros(67))
plt.scatter(data.loc[data['bio_rec'] == 1]['noHypNec_vascDensNoPref_038_ADCT2w'].to_numpy(),
            np.ones(9))

ax.set(title = 'Density', xlabel = 'noHypNec_vascDensNoPref_05_ADC', ylabel = 'Bio_rec')
ax.legend(['No. bio. rec.', 'Bio. rec.'])

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
#logReg = RidgeClassifier()
#logReg = LinearSVC()

#tx = [data[['age', 'gleason', 'psa', 'max_tum_area', 'tum_vol', 'total_dose',
#            'ADC_ave', 'T2w_ave']],
#    data[['max_tum_area', 'total_dose', 'ADC_ave']],
#    data[['max_tum_area', 'total_dose', 'ADC_ave', 'vascDensUniPref_03_ADC']],
#    data[['max_tum_area', 'total_dose', 'ADC_ave', 'vascDensUniPref_03_simp_ADC']],
#    data[['max_tum_area', 'total_dose', 'ADC_ave', 'vascDensUniPref_038_ADC']],
#    data[['max_tum_area', 'total_dose', 'ADC_ave', 'vascDensUniPref_038_simp_ADC']],
#    data[['max_tum_area', 'total_dose', 'ADC_ave', 'vascDensNoPref_05_ADC']],
#    data[['max_tum_area', 'total_dose', 'ADC_ave', 'vascDensNoPref_05_simp_ADC']],
#    data[['max_tum_area', 'total_dose', 'ADC_ave', 'noHypNec_vascDensNoPref_03_ADC']],
#    data[['max_tum_area', 'total_dose', 'ADC_ave', 'noHypNec_vascDensNoPref_03_simp_ADC']],
#    data[['max_tum_area', 'total_dose', 'ADC_ave', 'noHypNec_vascDensNoPref_05_ADC']],
#    data[['max_tum_area', 'total_dose', 'ADC_ave', 'noHypNec_vascDensNoPref_05_simp_ADC']]]
    
#tx = [data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave']],
#      data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave', 'vascDensUniPref_038_ADCT2w']]]
    
#tx = [data[['max_tum_area', 'total_dose', 'ADC_ave']],
#      data[['vascDensUniPref_038_ADC']],
#      data[['vascDensUniPref_038_simp_ADC']],
#      data[['vascDensUniPref_03_ADC']],
#      data[['vascDensUniPref_03_simp_ADC']],
#      data[['lowAlphaG1G2_vascDensUniPref_038_simp_ADC']]]

#tx = [data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave']],
#      data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave', 'vascDensUniPref_038_ADCT2w']],
#      data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave', 'vascDensUniPref_038_simp_ADCT2w']],
#      data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave', 'vascDensUniPref_05_ADCT2w']],
#      data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave', 'lowTTum_vascDensUniPref_038_ADCT2w']],
#      data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave', 'lowAlphaG1G2_vascDensUniPref_038_ADCT2w']],
#      data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave', 'lowPO2Nec_vascDensUniPref_038_ADCT2w']]]

#tx = [data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave']],
#      data[['vascDensNoPref_05_ADC']],
#      data[['vascDensNoPref_05_simp_ADC']],
#      data[['noHypNec_vascDensNoPref_05_ADC']],
#      data[['noHypNec_vascDensNoPref_05_simp_ADC']]]

tx = [data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave']],
      data[['init_tum_area']],
      data[['noHypNec_vascDensNoPref_038_ADCT2w']],
      data[['noHypNec_vascDensNoPref_038_simp_ADCT2w']]]

#tx = [data[['max_tum_area']],
#      data[['total_dose']],
#      data[['ADC_ave']],
#      data[['T2w_ave']],
#      data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave']],
#      data[['vascDensUniPref_038_ADCT2w']],
#      data[['max_tum_area', 'total_dose', 'ADC_ave']],
#      data[['vascDensUniPref_038_ADC']],
#      data[['max_tum_area', 'total_dose', 'T2w_ave']]]
    
y = data[['bio_rec']].to_numpy().ravel() 


#oversample = SMOTE()
#xO, yO = oversample.fit_resample(tx[0], y)

scores = np.zeros((N, K, len(tx)))
    
for j in range(N) :
        cv = StratifiedKFold(n_splits = K, shuffle = True)
        for i, x in enumerate(tx) :
            scores[j, :, i] = cross_val_score(logReg, x, y, cv = cv,
                  scoring = 'roc_auc')
        
scoresMean = np.mean(scores, axis = 1)
#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['all_features',
#                                     '3_features',
#                                     'vascDensUniPref_03_ADC',
#                                     'vascDensUniPref_03_simp_ADC',
#                                     'vascDensUniPref_038_ADC',
#                                     'vascDensUniPref_038_simp_ADC',
#                                     'vascDensNoPref_05_ADC',
#                                     'vascDensNoPref_05_simp_ADC',
#                                     'noHypNec_vascDensNoPref_03_ADC',
#                                     'noHypNec_vascDensNoPref_03_simp_ADC',
#                                     'noHypNec_vascDensNoPref_05_ADC',
#                                     'noHypNec_vascDensNoPref_05_simp_ADC'])
    
#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['3_features_ADC',
#                                     'vascDensUniPref_038_ADC',
#                                     'vascDensUniPref_038_simp_ADC',
#                                     'vascDensUniPref_03_ADC',
#                                     'vascDensUniPref_03_simp_ADC',
#                                     'lowAlphaG1G2_vascDensUniPref_038_simp_ADC'])
    
scoresMean = pd.DataFrame(data = scoresMean,
                          columns = ['4_features',
                                     'init_tum_area',
                                     'noHypNec_vascDensNoPref_038_ADCT2w',
                                     'noHypNec_vascDensNoPref_038_simp_ADCT2w'])

#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['4_features',
#                                     'vascDensUniPref_038_ADCT2w',
#                                     'vascDensUniPref_038_simp_ADCT2w',
#                                     'vascDensUniPref_05_ADCT2w',
#                                     'lowTTum_vascDensUniPref_038_ADCT2w',
#                                     'lowAlphaG1G2_vascDensUniPref_038_ADCT2w',
#                                     'lowPO2Nec_vascDensUniPref_038_ADCT2w'])

#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['4_features',
#                                    'vascDensUniPref_038_ADCT2w'])


#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['4_features',
#                                     'vascDensNoPref_05_ADC',
#                                     'vascDensNoPref_05_simp_ADC',
#                                     'noHypNec_vascDensNoPref_05_ADC',
#                                     'noHypNec_vascDensNoPref_05_simp_ADC'])
    
#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['max_tum_area',
#                                     'total_dose',
#                                     'ADC_ave',
#                                     'T2w_ave',
#                                     '4_features',
#                                     'vascDensUniPref_038_ADCT2w',
#                                     '3_features_ADC',
#                                     'vascDensUniPref_038_ADC',
#                                     '3_features_T2w'])

print('Mean AUC')
print(scoresMean.mean(axis = 0))

print('Standard deviation AUC')
print(scoresMean.std(axis = 0))

print('Median AUC')
print(scoresMean.median(axis = 0))
    
#%%
plt.close('all')
tcolor = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
ax = sns.boxplot(data = scoresMean, orient = 'v', palette = tcolor)
#ax = sns.swarmplot(data = scoresMean, color = ".25")
#add_stat_annotation(ax, data = scoresMean,
#                    box_pairs = [('vascDensUniPref_038_ADC',
#                                 'vascDensUniPref_038_simp_ADC'),
#                                 ('vascDensUniPref_03_ADC',
#                                 'vascDensUniPref_03_simp_ADC')],
#                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
#                    line_offset = -0.1, verbose = 2)
#add_stat_annotation(ax, data = scoresMean,
#                    box_pairs = [('vascDensUniPref_038_ADC', 
#                                'vascDensUniPref_03_ADC'),
#                                ('vascDensUniPref_038_ADC', 
#                                'lowAlphaG1G2_vascDensUniPref_038_simp_ADC')],
#                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
#                    line_offset = 0.0, verbose = 2)
#                                
#add_stat_annotation(ax, data = scoresMean,
#                    box_pairs = [('3_features_ADC', 
#                                'vascDensUniPref_038_ADC'),
#                                 ('3_features_ADC', 
#                                'vascDensUniPref_03_ADC'),
#                                ('3_features_ADC', 
#                                'lowAlphaG1G2_vascDensUniPref_038_simp_ADC')],
#                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
#                    line_offset = -0.05, verbose = 2)

add_stat_annotation(ax, data = scoresMean,
                    box_pairs = [('4_features',
                                'noHypNec_vascDensNoPref_038_ADCT2w'),
                                ('noHypNec_vascDensNoPref_038_ADCT2w',
                                'noHypNec_vascDensNoPref_038_simp_ADCT2w')],
                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
                    line_offset = -0.05, verbose = 2)

#add_stat_annotation(ax, data = scoresMean,
#                    box_pairs = [('4_features', 
#                                'vascDensUniPref_038_ADCT2w'),
#                                ('vascDensUniPref_038_ADCT2w', 
#                                'vascDensUniPref_038_simp_ADCT2w'),
#                                ('4_features',
#                                 'lowTTum_vascDensUniPref_038_ADCT2w'),
#                                ('4_features',
#                                 'lowAlphaG1G2_vascDensUniPref_038_ADCT2w'),
#                                ('4_features',
#                                 'lowPO2Nec_vascDensUniPref_038_ADCT2w')],
#                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
#                    line_offset = -0.05, verbose = 2)
#                                
#add_stat_annotation(ax, data = scoresMean,
#                    box_pairs = [('4_features', 
#                                'vascDensNoPref_05_ADC'),
#                                ('4_features', 
#                                'noHypNec_vascDensNoPref_05_ADC'),
#                                ('vascDensNoPref_05_ADC', 
#                                'vascDensNoPref_05_simp_ADC'),
#                                ('noHypNec_vascDensNoPref_05_ADC', 
#                                'noHypNec_vascDensNoPref_05_simp_ADC')],
#                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
#                    line_offset = -0.05, verbose = 2)
                               
#add_stat_annotation(ax, data = scoresMean,
#                    box_pairs = [('4_features', 
#                                'vascDensUniPref_038_ADCT2w'),
#                                ('4_features', 
#                                '3_features_ADC'),
#                                ('4_features', 
#                                '3_features_T2w'),
#                                ('3_features_ADC', 
#                                'vascDensUniPref_038_ADC')],
#                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
#                    line_offset = -0.05, verbose = 2)
                                
ax.set(ylabel = 'AUC', ylim = [0.5, 1])
#ax.set_xticklabels(['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave',
#                    '4_feat.', 'ref_ADCT2w', '3_feat._ADC', 'ref_ADC', 
#                    '3_feat._T2w'], rotation = 45, ha = 'center')

ax.set_xticklabels([], ha = 'center')

#%%   
plt.close('all')
plt.rcParams.update({'font.size': 32})  
N = 1000
K = 3

#tx = [data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave']],
#      data[['noHypNec_vascDensNoPref_038_ADCT2w']],
#      data[['noHypNec_vascDensNoPref_038_simp_ADCT2w']]]

#tx = [data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave']],
#      data[['init_tum_area']],
#      data[['noHypNec_vascDensNoPref_038_ADCT2w']],
#      data[['noHypNec_vascDensNoPref_038_simp_ADCT2w']]]

tx = [data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave']],
      data[['init_tum_area_ADC']],
      data[['init_tum_area_ADCT2w']]]

#
#tx = [data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave']],
#      data[['vascDensNoPref_05_ADC']]]

#tx = [data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave']],
#      data[['vascDensUniPref_038_ADCT2w']],
#      data[['vascDensUniPref_038_ADCT2w']],
#      data[['vascDensUniPref_038_simp_ADCT2w']],
#      data[['lowTTum_vascDensUniPref_038_ADCT2w']],
#      data[['lowAlphaG1G2_vascDensUniPref_038_ADCT2w']],
#      data[['lowPO2Nec_vascDensUniPref_038_ADCT2w']]]

#tx = [data[['max_tum_area', 'total_dose', 'ADC_ave']],
#      data[['vascDensUniPref_038_ADC']],
#      data[['vascDensUniPref_038_simp_ADC']],
#      data[['vascDensUniPref_03_ADC']],
#      data[['vascDensUniPref_03_simp_ADC']],
#      data[['lowAlphaG1G2_vascDensUniPref_038_simp_ADC']],
#      data[['noHypNec_vascDensNoPref_05_ADC']],
#      data[['noHypNec_vascDensNoPref_05_simp_ADC']]]

#tx = [data[['max_tum_area', 'total_dose', 'ADC_ave']],
#      data[['vascDensUniPref_03_ADC']],
#      data[['vascDensUniPref_03_simp_ADC']],
#      data[['vascDensUniPref_038_ADC']],
#      data[['vascDensUniPref_038_simp_ADC']],
#      data[['noHypNec_vascDensNoPref_03_ADC']],
#      data[['noHypNec_vascDensNoPref_03_simp_ADC']]]
y = data[['bio_rec']].to_numpy().ravel()  

mean_fpr = np.linspace(0, 1, 500)
mean_tprK = np.zeros((N, 500, len(tx)))
auc_ = np.zeros((N, len(tx))) 

logReg = LogisticRegression()
#logReg = LinearSVC()

fig, ax = plt.subplots() 
tcolor = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
for j in range(N) :
    cv = StratifiedKFold(n_splits = K, shuffle = True)    
    for i, x in enumerate(tx) :
        x = x.to_numpy()
        for k, (train, test) in enumerate(cv.split(x, y)) :
            logReg.fit(x[train], y[train])
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
for i in range(len(tx)) :
    ax.plot(mean_fpr, mean_tprN[:, i], linewidth = 6, color = tcolor[i])
    ci = stats.t.interval(level, dof, mean_tprN[:, i], std_tprN[:, i])
    ax.fill_between(mean_fpr, ci[0], ci[1], color = tcolor[i],
                    alpha = 0.1)

ax.plot([0, 1], [0, 1], '--', color = 'tab:gray', linewidth = 6)
ax.set(xlabel = 'FPR', ylabel = 'TPR', ylim = [0, 1])
mean_auc = np.mean(auc_, axis = 0)
std_auc = np.std(auc_, axis = 0)
ax.legend(["Total dose and pre-treat.\nimaging parameters\n"
           "(AUC = %0.2f $\pm$ %0.2f)" % (mean_auc[0], std_auc[0]),
           "Tum. area at 8 w. from\ncomprehensive model\n"
           "(AUC = %0.2f $\pm$ %0.2f)" % (mean_auc[1], std_auc[1]),
           "Tum. area at 8 w. from\nreduced model\n"
           "(AUC = %0.2f $\pm$ %0.2f)" % (mean_auc[2], std_auc[2])],
#           '(AUC = %0.2f)' % mean_auc[3],
#           '(AUC = %0.2f)' % mean_auc[4],
#           '(AUC = %0.2f)' % mean_auc[5],
#           '(AUC = %0.2f)' % mean_auc[6]],
           loc = 'lower right')
#
#           "Tum. area at 0 w.\n"
#           "(AUC = %0.2f $\pm$ %0.2f)" % (mean_auc[1], std_auc[1]),
#ax.legend(["Total dose and pre-treat.\nimaging parameters\n"
#           "(AUC = %0.3f $\pm$ %0.3f)" % (mean_auc[0], std_auc[0]),
#           "Initial tum. area\n(AUC = %0.3f $\pm$ %0.3f)" % (mean_auc[1],
#                                std_auc[1]),
#           "Tum. area at 8 w. from\ncomprehensive model\n"
#           "(AUC = %0.3f $\pm$ %0.3f)" % (mean_auc[2], std_auc[2]),
#           "Tum. area at 8 w. from\nreduced model\n"
#           "(AUC = %0.3f $\pm$ %0.3f)" % (mean_auc[3], std_auc[3])],
##           '(AUC = %0.2f)' % mean_auc[3],
##           '(AUC = %0.2f)' % mean_auc[4],
##           '(AUC = %0.2f)' % mean_auc[5],
##           '(AUC = %0.2f)' % mean_auc[6]],
#           loc = 'lower right')
#
_, p = stats.wilcoxon(auc_[:, 0], auc_[:, 1], alternative = 'less')
print("4_features, comp.:", p)
_, p = stats.wilcoxon(auc_[:, 0], auc_[:, 2], alternative = 'less')
print("4_features, red.:", p)
_, p = stats.wilcoxon(auc_[:, 1], auc_[:, 2])
print("comp., red.:", p)

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

tx = [data[['age', 'gleason', 'psa', 'max_tum_area', 'tum_vol', 'total_dose',
            'ADC_ave', 'T2w_ave']],
      data[['max_tum_area', 'total_dose', 'ADC_ave']],
      data[['max_tum_area', 'total_dose', 'ADC_ave', 'T2w_ave']],
      data[['two_mon_tum_area']]]
y = data[['bio_rec']].to_numpy().ravel() 

scores = np.zeros((N, K, len(tx)))
    
for j in range(N) :
        cv = StratifiedKFold(n_splits = K, shuffle = True)
        for i, x in enumerate(tx) :
            scores[j, :, i] = cross_val_score(logReg, x, y, cv = cv,
                  scoring = 'roc_auc')
        
scoresMean = np.mean(scores, axis = 1)
scoresMean = pd.DataFrame(data = scoresMean,
                          columns = ['all_features',
                                     '3_features',
                                     '4_features',
                                     'two_mon_tum_area'])
#%%
plt.close('all')
ax = sns.boxplot(data = scoresMean, orient = 'v', palette = 'Set2')
add_stat_annotation(ax, data = scoresMean,
                    box_pairs=[('all_features', 'two_mon_tum_area'),
                               ('3_features', 'two_mon_tum_area'),
                               ('3_features', '4_features')],
                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
                    line_offset = -0.05, verbose = 2)
#_, p = stats.wilcoxon(scoresMean[['all_features']].to_numpy().ravel(),
#                      scoresMean[['3_features']].to_numpy().ravel())
#print(p)
ax.set(title = 'AUC', ylabel = 'AUC', ylim = [None, 1])
ax.set_xticklabels(['All features',
                    'Initial ADC and\nmaximal tumor area\nand total dose',
                    '4',
                    '$In$ $silico$ 2 months\ntumor area'], ha = 'center')

#%%  
plt.close('all')
plt.rcParams.update({'font.size': 32})  
N = 1000
K = 3
#data[['two_mon_tum_area', 'three_mon_tum_area', 'five_mon_tum_area',
#      'int_two_mon_tum_area', 'int_five_mon_tum_area']].to_numpy(),
tx = [data[['age', 'gleason', 'max_tum_area', 'tum_vol', 'total_dose',
            'ADC_ave', 'T2w_ave']].to_numpy(),
      data[['ADC_ave', 'max_tum_area', 'total_dose']].to_numpy(),
      data[['two_mon_tum_area']].to_numpy(),
      data[['age', 'gleason', 'max_tum_area', 'tum_vol', 'total_dose',
            'ADC_ave', 'T2w_ave', 'two_mon_tum_area']].to_numpy(),]
y = data[['bio_rec']].to_numpy().ravel()  

mean_fpr = np.linspace(0, 1, 100)
mean_tprK = np.zeros((N, 100, len(tx)))
auc_ = np.zeros((N, len(tx))) 

logReg = LogisticRegression()

fig, ax = plt.subplots() 

for j in range(N) :
    cv = StratifiedKFold(n_splits = K, shuffle = True)    
    for i, x in enumerate(tx) :
        for k, (train, test) in enumerate(cv.split(x, y)) :
            logReg.fit(x[train], y[train])
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
ax.plot(mean_fpr, mean_tprN, linewidth = 6)
tcolor = ['tab:blue', 'tab:orange', 'tab:green', 'tab:purple']
for i in range(len(tx)) :
    ci = stats.t.interval(level, dof, mean_tprN[:, i], std_tprN[:, i])
    ax.fill_between(mean_fpr, ci[0], ci[1], color = tcolor[i],
                    alpha = 0.1)

ax.plot([0, 1], [0, 1], '--', color = 'tab:gray', linewidth = 6)
ax.set(title = 'ROC', xlabel = 'FPR', ylabel = 'TPR', ylim = [0, 1])
mean_auc = np.mean(auc_, axis = 0)
ax.legend(['All features (AUC = %0.2f)' % mean_auc[0],
           'Initial ADC and maximal\ntumor area and total dose\n(AUC = %0.2f)'
           % mean_auc[1],
           '$In$ $silico$ 2 months tumor\narea (AUC = %0.2f)' % mean_auc[2],
           'All features and $in$ $silico$\n2 months tumor area\n(AUC = %0.2f)'
           % mean_auc[3]],
           loc = 'lower right')

_, p = stats.wilcoxon(auc_[:, 1], auc_[:, 0])
print(p)
_, p = stats.wilcoxon(auc_[:, 0], auc_[:, 2], alternative = 'less')
print(p)
_, p = stats.wilcoxon(auc_[:, 1], auc_[:, 2], alternative = 'less')
print(p)
_, p = stats.wilcoxon(auc_[:, 0], auc_[:, 3], alternative = 'less')
print(p)


#Con un hypNecThres de 0.7 y una densidad vascular de 0.05, en un tumor de tumDens Area =1
#no aparecen regiones hipoxicas. Es como si el umbral estuviera a 0
#Si se baja la densidad vascular a 0.03, aparecen algunas regiones hipoxicas,
# aunaue no necesariamente en el centro del tumor, ya que la distribucion de las
#celulas endoteliales es uniforme. Con todo, parece la opcion mas realista.
#Lanzar una simulacion con densVasc = 0.03 (teniendo prioridad las celulas
#endoteliales y hypNecThres = 0.7)
#Repetir con estos valores para el modelo simplificado  

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

adcRec = data.loc[data['bio_rec'] == 1][['noHypNec_vascDensNoPref_038_simp_ADCT2w']].to_numpy()
adcRec = adcRec.ravel()
ecdfRec = ECDF(adcRec)
popt, pcov = curve_fit(fsigmoid, ecdfRec.x[1:], ecdfRec.y[1:], p0 = [1, 50])
for j in range(301):
    sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
ax.plot(ecdfRec.x, ecdfRec.y, linewidth = 6)
#ax.plot(x, sigmoid, linewidth = 6)

adcNoRec = data.loc[data['bio_rec'] == 0][['noHypNec_vascDensNoPref_038_simp_ADCT2w']].to_numpy()
adcNoRec = adcNoRec.ravel()
ecdfNoRec = ECDF(adcNoRec)
popt, pcov = curve_fit(fsigmoid, ecdfNoRec.x[1:], ecdfNoRec.y[1:],
                       p0 = [1, 50])
for j in range(301):
    sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
ax.plot(ecdfNoRec.x, ecdfNoRec.y, linewidth = 6)
#ax.plot(x, sigmoid, linewidth = 6)
ax.set(xlabel = 'average ADC', ylabel = 'cumulative density', title = 'ADC')
ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])

thresholds = np.linspace(float(data[['noHypNec_vascDensNoPref_038_ADCT2w']].min()),
                         float(data[['noHypNec_vascDensNoPref_038_ADCT2w']].max()), 20)

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
plt.close('all')
#tx = ['ADC_ave', 'max_tum_area', 'tum_vol', 'T2w_ave',
#      'epi_dens_ADC_ave', 'epi_dens_T2w_ave', 'epi_dens_mean']
#txNames = ['ADC', 'maxTumArea', 'tumVol', 'T2w', 'dens(ADC)', 'dens(T2w)',
#           'dens(ADC&T2w)']

#tx = ['ADC_ave', 'ADC_ave_norm', 'epi_dens_norm_ADC_ave',
#      'epi_dens_ADC_ave', 'T2w_ave', 'T2w_average_norm',
#      'epi_dens_norm_T2w_ave', 'epi_dens_T2w_ave', 'epi_dens_mean',
#      '5_mon_tum_area']
#txNames = ['ADC', 'ADCNorm', 'densNorm(ADC)', 'dens(ADC)', 'T2w', 'T2wNorm', 
#           'densNorm(T2w)', 'dens(T2w)', 'dens(ADC&T2w)', '5MonTumArea']

tx = ['ADC_ave', 'epi_dens_ADC_ave', 'max_tum_area', 'five_mon_tum_area']
txNames = ['ADC', 'dens(ADC)', 'maxTumArea','5MonTumArea']


N = 1000
K = 5
scores = np.zeros((N, K, len(tx)))
logReg = LogisticRegression()   
y = data[['bio_rec']].to_numpy().ravel() 
    
for i, el in enumerate(tx) :    
    x = data[[el]].to_numpy()
    
    fig, ax = plt.subplots() 
    mean_fpr = np.linspace(0, 1, 100)
    mean_tprK = np.zeros((N, 100))   
    logReg = LogisticRegression()
    
    for j in range(N) :
        cv = StratifiedKFold(n_splits = K, shuffle = True)
                 
#        fig, ax = plt.subplots() 
        for k, (train, test) in enumerate(cv.split(x, y)) :
            logReg.fit(x[train], y[train])
            probas = logReg.predict_proba(x[test])
            fpr, tpr, _= roc_curve(y[test], probas[:, 1])   
            mean_tprK[j, :] += np.interp(mean_fpr, fpr, tpr)
            mean_tprK[j, 0] = 0.0

#            ax.plot(fprK, tprK, linewidth = 6)
            
        mean_tprK[j, :] /= K
        mean_tprK[j, -1] = 1.0
        
#        ax.plot(mean_fpr, mean_tpr[j, :], linewidth = 6)
    
    print(logReg.classes_)
    mean_tprN = np.mean(mean_tprK, axis = 0)
    std_tprN = np.std(mean_tprK, axis = 0)
    auc_ = np.trapz(mean_tprN, mean_fpr)
    level = 0.95
    dof = K - 1

    ci = stats.t.interval(level, dof, mean_tprN, std_tprN)

    ax.plot(mean_fpr, mean_tprN, linewidth = 6)
    ax.fill_between(mean_fpr, ci[0], ci[1], color = 'tab:blue', alpha = 0.1)
    ax.plot([0, 1], [0, 1], '--', color = 'tab:gray', linewidth = 6)
    ax.set(title = 'ROC - ' + txNames[i], xlabel = 'FPR', ylabel = 'TPR',
           ylim = [0, 1])
    ax.legend(['ROC (AUC = %0.2f)' % auc_], loc = 'lower right')

#%%
plt.close('all')
#tx = ['ADC_average', 'max_tumor_area', 'tumor_volume', 'T2w_average',
#      'epi_dens_ADC_ave', 'epi_dens_T2w_ave', 'epi_dens_mean']
#txNames = ['ADC', 'maxTumArea', 'tumVol', 'T2w', 'dens(ADC)', 'dens(T2w)',
#           'dens(ADC&T2w)']

#tx = ['ADC_average', 'max_tumor_area', 'T2w_average']
#txNames = ['average ADC', 'max_tumor_area', 'average T2w']

#tx = ['ADC_average', 'ADC_average_norm', 'epi_dens_norm_ADC_ave',
#      'epi_dens_ADC_ave', 'T2w_average', 'T2w_average_norm',
#      'epi_dens_norm_T2w_ave', 'epi_dens_T2w_ave', 'epi_dens_mean', 'random']
#txNames = ['ADC', 'ADCNorm', 'densNorm(ADC)', 'dens(ADC)', 'T2w', 'T2wNorm', 
#           'densNorm(T2w)', 'dens(T2w)', 'dens(ADC&T2w)', 'random']

tx = ['epi_dens_ADC_ave', 'epi_dens_T2w_ave', 'epi_dens_mean', 'random']
txNames = ['dens(ADC)', 'dens(T2w)', 'dens(ADC&T2w)', 'random']

N = 1000
K = 3
scores = np.zeros((N, K, len(tx)))
logReg = LogisticRegression()
    
for i, el in enumerate(tx) :    
    x = data[[el]].to_numpy()
    y = data[['bio_rec']].to_numpy().ravel() 
    
    
    
    fig, ax = plt.subplots() 
    mean_fpr = np.linspace(0, 1, 100)
    mean_tprK = np.zeros((N, 100))   
    
    for j in range(N) :
        cv = StratifiedKFold(n_splits = K, shuffle = True)
                 
#        fig, ax = plt.subplots() 
        for k, (train, test) in enumerate(cv.split(x, y)) :
            logReg.fit(x[train], y[train])
            probas = logReg.predict_proba(x[test])
            fpr, tpr, _= roc_curve(y[test], probas[:, 0],
                                 drop_intermediate = False) 
            
            mean_tprK[j, :] += np.interp(mean_fpr, fpr, tpr)
            mean_tprK[j, 0] = 0.0

#            ax.plot(fprK, tprK, linewidth = 6)
            
        mean_tprK[j, :] /= K
        mean_tprK[j, -1] = 1.0
        
#        ax.plot(mean_fpr, mean_tpr[j, :], linewidth = 6)
    
    mean_tprN = np.mean(mean_tprK, axis = 0)
    std_tprN = np.std(mean_tprK, axis = 0)
    auc_ = np.trapz(mean_tprN, mean_fpr)
    level = 0.95
    dof = K - 1

    ci = stats.t.interval(level, dof, mean_tprN, std_tprN)

    ax.plot(mean_fpr, mean_tprN)
    ax.fill_between(mean_fpr, ci[0], ci[1], color = 'tab:blue', alpha = 0.1)
    ax.plot([0, 1], [0, 1], '--', color = 'tab:gray', linewidth = 6)
    ax.set(title = 'ROC - ' + txNames[i], xlabel = 'FPR', ylabel = 'TPR')
    ax.legend(['ROC (AUC = %0.2f)' % auc_], loc = 'lower right')

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
x = data[['ADC_average', 'T2w_average']].to_numpy()
y = data[['bio_rec']].to_numpy().ravel()
xtrain, xtest, ytrain, ytest = train_test_split(x, y)

logReg = LogisticRegression()
logReg.fit(xtrain, ytrain)
ypred = logReg.predict(xtest)

confusion_matrix(ytest, ypred)

probas = logReg.predict_proba(xtest)
fpr, tpr, thresholds = roc_curve(ytest, probas[:, 0],
                                 pos_label = logReg.classes_[0],
                                 drop_intermediate = False)
rocAuc = auc(fpr, tpr)

fig, ax = plt.subplots()
ax.plot(thresholds, fpr, linewidth = 6)
ax.plot(thresholds, tpr, linewidth = 6)
ax.set(title = 'Rates - ' + txNames[i], xlabel = 'threshold',
       ylabel = 'rate')
ax.legend(['False positives', 'True positives'])

fig, ax = plt.subplots()    
ax.plot(fpr, tpr, linewidth = 6)
ax.plot([0, 1], [0, 1], '--k', linewidth = 6)
ax.set(title = 'ROC - average ADC and average T2w')
ax.legend(['ROC (AUC = %0.2f)' % rocAuc])

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
X = np.loadtxt('Tumor_76_ADC_NoHormono_NoDates.dat')
nSamples = np.size(X, 0)
y = np.loadtxt('recidive.dat')
plt.close('all')
nfig = 0

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
B = 10
nNeigh = 10
RMSE = np.zeros([nNeigh, B])
for b in range(nNeigh):
    itrain, itest = train_test_split(range(nSamples - 1), test_size = 0.2)
    Xtrain = X[itrain, :]
    ytrain = y[np.asarray(itrain)] 
    Xtest = X[itest, :]
    ytest = y[np.asarray(itest)]
    for i in range(1, nNeigh + 1):
        knn = KNeighborsClassifier(n_neighbors = nNeigh, weights = 'distance');
        knn.fit(Xtrain, ytrain)
        ypred = knn.predict(Xtest)
        RMSE[i - 1, b] = np.sqrt(np.mean((ytest - ypred)**2))
nfig += 1
plt.figure(nfig)
plt.boxplot(RMSE.T, labels = ['1', '2', '3', '4', '5', '6', '7', '8', '9', 
                              '10'])
plt.title("kNN")
plt.ylim(0.0, 0.7)
plt.ylabel('RMSE')
plt.xlabel('Number of neighbors')

knn = KNeighborsClassifier(n_neighbors = 2, weights = 'distance')
knn.fit(Xtrain, ytrain)
ypred = knn.predict(Xtest)

nfig += 1
plt.figure(nfig)
plt.subplot(2, 1, 1)
plt.title("kNN")
line1, = plt.plot(range(len(ytest)), ytest,'-k') 
line2, = plt.plot(range(len(ytest)), ypred,'*-r') 
plt.ylim([-2.5, 2.5])
plt.legend([line1, line2], ['True ENSO index', 'Estimated ENSO index'])
plt.subplot(2, 1, 2)
plt.plot(ypred, ytest,'*k')
plt.plot([-2, 2],[-2, 2],'-r')
plt.xlim([-2.5, 2.5])
plt.ylim([-2.5, 2.5])
plt.xlabel('True ENSO index')
plt.ylabel('Estimated ENSO index')

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