#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.ensemble import RandomSurvivalForest
from sklearn.model_selection import (train_test_split, StratifiedKFold,
                                     LeaveOneOut, StratifiedShuffleSplit)
from sksurv.metrics import (concordance_index_censored, concordance_index_ipcw,
                            cumulative_dynamic_auc)
import seaborn as sns
from statannot import add_stat_annotation
    
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
#data = pd.read_csv('../../Carlos/Results/Recurrence/noHypNec_vascDensNoPref0.05/'
#                   'noHypNec_vascDensNoPref_05.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/vascDensUniPref0.03_ADC/'
#                   'vascDensUniPref_03.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/rec_summary_8wTumVol.csv')
data = pd.read_csv('../../Carlos/Results/Recurrence/simp/rec_summary_8wTumVol_norm.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/simp/rec_summary_8wIntTumVol.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120/TTum330_alphaG1120.csv')

#%%
cph = CoxPHFitter(penalizer = 0.1, l1_ratio = 0.0)   
cph.fit(data, duration_col = 'bio_rec_6_delay', event_col = 'bio_rec_6')
cph.print_summary() 
cph.plot()  
#cph.predict_survival_function(data)
#cph.plot_partial_effects_on_outcome(covariates = 'gleason', values = [6, 7, 8],
#                                    cmap='coolwarm')


#%%
N = 100
#x = data
#x = x.drop(columns = ['bio_rec_6', 'bio_rec_6_delay'])
tx = [data[['ADC_med', 'T2w_diff_var_mean', 'tum_vol', 'T2w_contrast_mean']],
      data[['tum_vol', 'ADC_ave', 'T2w_ave']],
      data[['tum_area_from_vol', 'dens_ADCT2w']],
      data[['init_tum_area']],
      data[['TTum330_alphaG1120']]]

y = np.zeros(76, dtype={'names':('bio_rec_6', 'bio_rec_6_delay'),
                          'formats':('bool', 'int')})
y['bio_rec_6'] = data[['bio_rec_6']].to_numpy().ravel()
y['bio_rec_6_delay'] = data[['bio_rec_6_delay']].to_numpy().ravel()

nTest = 38
auc = np.zeros((len(tx), 3))
meanAuc = np.zeros((N, len(tx)))
cph = CoxPHSurvivalAnalysis()
sss = StratifiedShuffleSplit(n_splits = N, test_size = nTest)

for i, x in enumerate(tx) :
    x = x.to_numpy()
    j = 0
    for train, test in sss.split(x, y['bio_rec_6']):
        xtrain, ytrain = x[train], y[train]
        xtest, ytest = x[test], y[test]
        cph.fit(xtrain, ytrain)
        ypred = cph.predict(xtrain)
        time = np.linspace(np.min(ytrain['bio_rec_6_delay']) + 1,
                           np.max(ytrain['bio_rec_6_delay']) - 1, 3)
        auc[i, :], meanAuc[j, i] = cumulative_dynamic_auc(ytrain, ytrain,
           ypred, time)
        j += 1

meanAuc = pd.DataFrame(data = meanAuc, 
                       columns = ['4_feat.', "Prediction 1", 'Prediction 2',
                                  'Prediction 3', 'Prediction 4'])
meanMeanAuc = meanAuc.mean(axis = 0)
stdMeanAuc = meanAuc.std(axis = 0)

#%%
N = 100
#x = data
#x = x.drop(columns = ['bio_rec_6', 'bio_rec_6_delay'])
tx = [data[['ADC_med', 'T2w_diff_var_mean', 'tum_vol', 'T2w_contrast_mean']],
      data[['tum_vol', 'ADC_ave', 'T2w_ave']],
      data[['tum_area_from_vol', 'dens_ADCT2w']],
      data[['init_tum_area']],
      data[['TTum330_alphaG1120']]]

y = np.zeros(76, dtype={'names':('bio_rec_6', 'bio_rec_6_delay'),
                          'formats':('bool', 'int')})
y['bio_rec_6'] = data[['bio_rec_6']].to_numpy().ravel()
y['bio_rec_6_delay'] = data[['bio_rec_6_delay']].to_numpy().ravel()

time = np.arange(np.min(y['bio_rec_6_delay']) + 1,
                  np.max(y['bio_rec_6_delay']) - 1, 6)
auc = np.zeros((len(tx), len(time)))
meanAuc = np.zeros(len(tx))

cph = CoxPHSurvivalAnalysis()

for i, x in enumerate(tx) :
    x = x.to_numpy()
    cph.fit(x, y)
    ypred = cph.predict(x)
    auc[i, :], meanAuc[i] = cumulative_dynamic_auc(y, y, ypred, time)

#%%
plt.close('all')
fig, ax = plt.subplots() 
tcolor = [red, darkPurple, orangePurple, orange, greenBlue]

for i in range(len(tx)) :
    ax.plot(time, auc[i, :], marker = "o", markersize = 12, linewidth = 6,
            color = tcolor[i])

ax.set(xlabel = 't (months)', ylabel = 'Time-dependent AUC',
       ylim = [0.0, None])
ax.grid(True)
ax.legend(['Med. ADC, tum. vol., T2w diff. var.\nmean & '
           'T2w contrast mean\n(Mean AUC = %0.2f)' % meanAuc[0],
           'Prediction 1 (Mean AUC = %0.2f)' % meanAuc[1],
           'Prediction 2 (Mean AUC = %0.2f)' % meanAuc[2],
           'Prediction 3 (Mean AUC = %0.2f)' % meanAuc[3],
           'Prediction 4 (Mean AUC = %0.2f)' % meanAuc[4]],
          loc = 'lower right', fontsize = 24)

#%%
plt.close('all')
tcolor = [red, darkPurple, orangePurple, orange, greenBlue]
ax = sns.boxplot(data = meanAuc, orient = 'v', palette = tcolor)
#ax = sns.swarmplot(data = scores, color = ".25")

add_stat_annotation(ax, data = meanAuc,
                    box_pairs = [('4_feat.',
                                'Prediction 1'), 
                                ('Prediction 1',
                                'Prediction 2'),
                                ('Prediction 2',
                                'Prediction 3'),
                                ('Prediction 3',
                                'Prediction 4'),
                                ('4_feat.',
                                'Prediction 4')],
                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
                    line_offset = 0.025, verbose = 2)
                                
#add_stat_annotation(ax, data = scores,
#                    box_pairs = [('4_feat.',
#                                'Prediction 3'),
#                                ('Prediction 3',
#                                'vd02'),
#                                ('Prediction 3',
#                                'vd025'),
#                                ('Prediction 3',
#                                'vd03'),
#                                ('Prediction 3',
#                                'vd038'),
#                                ('Prediction 3',
#                                'vd05')],
#                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
#                    line_offset = 0.0, verbose = 2)

ax.set(ylabel = 'C-index', ylim = [.5, 1.])
ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_xticklabels(['4_feat.', "Prediction 1", 'Prediction 2',
                    'Prediction 3', 'Prediction 4'], ha = 'center') 
#ax.set_xticklabels(['4_feat.', "Prediction 3", 'vd = 2%', 'vd = 2.5%',
#                    'vd = 3%', 'vd = 3.8%', 'vd = 5%'], ha = 'center') 


#%%
kmf = KaplanMeierFitter()
T = data['bio_rec_6_delay']
E = data['bio_rec_6']
kmf.fit(T, event_observed = E)

fig, ax = plt.subplots() 
kmf.plot_survival_function(at_risk_counts = True)
plt.tight_layout()

#%%
#cov = 'age'
#nameCov = 'Age'
#unitCov = 'years'

#cov = 'tum_vol'
#nameCov = 'Tumor volume'
#unitCov = 'mm$^3$'

#cov = 'ADC_ave'
#nameCov = 'Average ADC'
#unitCov = ''

#cov = 'T2w_ave'
#nameCov = 'Average T2-w'
#unitCov = ''

#cov = 'T2w_diff_var_mean'
#nameCov = 'Mean T2-w difference variance'
#unitCov = ''

#cov = 'T2w_contrast_mean'
#nameCov = 'Mean T2-w contrast'
#unitCov = ''

#cov = 'tum_area_from_vol'
#nameCov = '$A_{tum}$'
#unitCov = ''

#cov = 'dens_ADCT2w'
#nameCov = r'$\rho_{tum}$'
#unitCov = ''

#cov = 'dens_ADC'
#nameCov = r'$\rho^{ADC}_{tum}$'
#unitCov = ''

#cov = 'init_tum_area'
#nameCov = '$\it{In}$ $\it{silico}$ initial n° of tumor cells'
#unitCov = ''

#cov = 'TTum330_alphaG1120'
#nameCov = '$\it{In}$ $\it{silico}$ n° of tumor cells at $t$ = 8 weeks'
#unitCov = ''

#cov = 'TTum330_alphaG1120_norm'
#nameCov = 'Norm. $\it{in}$ $\it{silico}$ n° of tumor cells at $t$ = 8 weeks'
#unitCov = ''

cov = 'TTum330_alphaG1120_probRF'
nameCov = 'Probability with RF from $\it{In}$ $\it{silico}$ n° of tumor cells at $t$ = 8 weeks'
unitCov = ''

thresholds = np.linspace(data[cov].min(), data[cov].max(), 100)
best_threshold = thresholds[0]
best_pvalue = 1.
for threshold in thresholds:
    grThres = data[cov] > threshold
    logRank = logrank_test(T[grThres], T[~grThres], E[grThres], E[~grThres],
                           alpha = .99)
    
    if (logRank.p_value < best_pvalue and len(T[grThres]) > 18 and
        len(T[~grThres]) > 18) :
        best_pvalue = logRank.p_value
        best_threshold = threshold


best_threshold = data[cov].median() 
grThres = data[cov] > best_threshold

#%%
fig, ax = plt.subplots() 
kmf.fit(T[grThres], event_observed = E[grThres], label = nameCov +
        ' $>$ %0.2f ' % best_threshold + unitCov + '\n(n = %i)' %
        len(T[grThres]))
#kmf.fit(T[grThres], event_observed = E[grThres], label = nameCov +
#        ' $>$ %0.0f ' % (best_threshold / (0.02 * 0.02)) + unitCov +
#        ' (n = %i)' % len(T[grThres]))
kmf.plot_survival_function(ax = ax, linewidth = 6)

kmf.fit(T[~grThres], event_observed = E[~grThres], label = nameCov +
        ' $\leq $ %0.2f ' % best_threshold + unitCov + '\n(n = %i)' %
        len(T[~grThres]))
#kmf.fit(T[~grThres], event_observed = E[~grThres], label = nameCov +
#        ' $\leq $ %0.0f ' % (best_threshold / (0.02 * 0.02)) + unitCov +
#        ' (n = %i)' % len(T[~grThres]))
kmf.plot_survival_function(ax = ax, linewidth = 6)

ax.set(xlabel = 't (months)',
       ylabel = 'Freedom from biochemical recurrence rate')
ax.legend(loc = 'lower left')

logRank = logrank_test(T[grThres], T[~grThres], E[grThres], E[~grThres],
                           alpha = .99)
logRank.print_summary()
ax.text(0.05, 0.4, 'p-value = %f' % logRank.p_value)

#%%
N = 100
K = 3

#tx = [data[['ADC_med', 'T2w_diff_var_mean', 'tum_vol', 'T2w_contrast_mean']],
#      data[['ADC_med', 'T2w_diff_var_mean', 'tum_vol', 'T2w_contrast_mean',
#            'ADC_ave', 'T2w_ave']],
#      data[['ADC_med', 'T2w_diff_var_mean', 'tum_vol', 'T2w_contrast_mean',
#            'tum_area_from_vol', 'dens_ADCT2w']],
#      data[['ADC_med', 'T2w_diff_var_mean', 'tum_vol', 'T2w_contrast_mean',
#            'init_tum_area']],
#      data[['ADC_med', 'T2w_diff_var_mean', 'tum_vol', 'T2w_contrast_mean',
#            'TTum330_alphaG1120']]]

tx = [data[['ADC_med', 'T2w_diff_var_mean', 'tum_vol', 'T2w_contrast_mean']],
      data[['tum_vol', 'ADC_ave', 'T2w_ave']],
      data[['tum_area_from_vol', 'dens_ADCT2w']],
      data[['init_tum_area']],
      data[['TTum330_alphaG1120_vd02']]]
      
y = np.zeros(76, dtype={'names':('bio_rec_6', 'bio_rec_6_delay'),
                          'formats':('bool', 'int')})
y['bio_rec_6'] = data[['bio_rec_6']].to_numpy().ravel()
y['bio_rec_6_delay'] = data[['bio_rec_6_delay']].to_numpy().ravel()

#rsf = RandomSurvivalForest(n_estimators = 100, min_samples_split = 10,
#                           min_samples_leaf = 15, max_features = "sqrt")
rsf = RandomSurvivalForest(n_estimators = 100, min_samples_leaf = 15)

scores = np.zeros((N, len(tx))) 

for j in range(N) :
    cv = StratifiedKFold(n_splits = K, shuffle = True)    
    for i, x in enumerate(tx) :
        x = x.to_numpy()
        for k, (train, test) in enumerate(cv.split(x, y['bio_rec_6'])) :
            xtrain, ytrain = x[train], y[train]
            xtest, ytest = x[test], y[test]
            rsf.fit(xtrain, ytrain)
            scores[j, i] += rsf.score(xtest, ytest)
scores /= K

mean_scores = np.mean(scores, axis = 0)
std_scores = np.std(scores, axis = 0, ddof = 1)

#level = 0.95
#dof = K - 1

scores = pd.DataFrame(data = scores,
                      columns = ['4_feat.', "Prediction 1", 'Prediction 2',
                                 'Prediction 3', 'Prediction 4'])

#scores = pd.DataFrame(data = scores,
#                      columns = ['4_feat.', 'Prediction 3', 'vd02', 'vd025',
#                                 'vd03', 'vd038', 'vd05'])
#%%
plt.close('all')
tcolor = [red, darkPurple, orangePurple, orange, greenBlue]
ax = sns.boxplot(data = scores, orient = 'v', palette = tcolor)
#ax = sns.swarmplot(data = scores, color = ".25")

add_stat_annotation(ax, data = scores,
                    box_pairs = [('4_feat.',
                                'Prediction 1'), 
                                ('Prediction 1',
                                'Prediction 2'),
                                ('Prediction 2',
                                'Prediction 3'),
                                ('Prediction 3',
                                'Prediction 4'),
                                ('4_feat.',
                                'Prediction 4')],
                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
                    line_offset = -0.05, verbose = 2)
                                
#add_stat_annotation(ax, data = scores,
#                    box_pairs = [('4_feat.',
#                                'Prediction 3'),
#                                ('Prediction 3',
#                                'vd02'),
#                                ('Prediction 3',
#                                'vd025'),
#                                ('Prediction 3',
#                                'vd03'),
#                                ('Prediction 3',
#                                'vd038'),
#                                ('Prediction 3',
#                                'vd05')],
#                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
#                    line_offset = 0.0, verbose = 2)

ax.set(ylabel = 'C-index', ylim = [.5, 1.])
ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_xticklabels(['4_feat.', "Prediction 1", 'Prediction 2',
                    'Prediction 3', 'Prediction 4'], ha = 'center') 
#ax.set_xticklabels(['4_feat.', "Prediction 3", 'vd = 2%', 'vd = 2.5%',
#                    'vd = 3%', 'vd = 3.8%', 'vd = 5%'], ha = 'center') 

    
#%%
N = 100
K = 3

x = data[['ADC_med', 'T2w_diff_var_mean', 'tum_vol', 'T2w_contrast_mean',
          'TTum330_alphaG1120']].to_numpy()
      
y = np.zeros(76, dtype={'names':('bio_rec_6', 'bio_rec_6_delay'),
                          'formats':('bool', 'int')})
y['bio_rec_6'] = data[['bio_rec_6']].to_numpy().ravel()
y['bio_rec_6_delay'] = data[['bio_rec_6_delay']].to_numpy().ravel()

#tNEstimators = [10, 20, 50, 100, 200, 500, 1000]
#tNSampleSplit = [2, 4, 6, 8, 10, 12, 14]
tNSampleLeaf = [3, 6, 9, 12, 15, 18, 21]

scores = np.zeros((N, len(tNSampleLeaf))) 

for j in range(N) :
    cv = StratifiedKFold(n_splits = K, shuffle = True)    
    for i, nSampleLeaf in enumerate(tNSampleLeaf) :
        rsf = RandomSurvivalForest(n_estimators = 100,
                                   min_samples_split = 10, 
                                   min_samples_leaf = nSampleLeaf,
                                   max_features = "sqrt")
        for k, (train, test) in enumerate(cv.split(x, y['bio_rec_6'])) :
            xtrain, ytrain = x[train], y[train]
            xtest, ytest = x[test], y[test]
            rsf.fit(xtrain, ytrain)
            scores[j, i] += rsf.score(xtest, ytest)
scores /= K

mean_scores = np.mean(scores, axis = 0)
std_scores = np.std(scores, axis = 0, ddof = 1)

#level = 0.95
#dof = K - 1

scores = pd.DataFrame(data = scores,
                      columns = ['10', '20', '50', '100', '200', '500',
                                 '1000'])
#%%
plt.close('all')
tcolor = [red, darkPurple, orangePurple, orange, greenBlue]
ax = sns.boxplot(data = scores, orient = 'v', palette = tcolor)
#ax = sns.swarmplot(data = scores, color = ".25")

add_stat_annotation(ax, data = scores,
                    box_pairs = [('10',
                                '20'), 
                                ('20',
                                '50'),
                                ('50',
                                '100'),
                                ('100',
                                '200'),
                                ('200',
                                '500'),
                                ('500',
                                '1000')],
                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
                    line_offset = -0.075, verbose = 2)
                                


ax.set(ylabel = 'C-index', ylim = [.5, 1.])
ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
#ax.set_xticklabels(['10', '20', '50', '100', '200', '500',
#                                 '1000'], ha = 'center') 

#ax.set_xticklabels([2, 4, 6, 8, 10, 12, 14], ha = 'center') 
ax.set_xticklabels([3, 6, 9, 12, 15, 18, 21], ha = 'center') 
