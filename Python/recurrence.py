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
#data = pd.read_csv('../../Carlos/Results/Recurrence/simp/rec_summary_8wTumVol.csv')
#data = pd.read_csv('../../Carlos/Results/Recurrence/simp/rec_summary_8wIntTumVol.csv')
data = pd.read_csv('../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120/TTum330_alphaG1120.csv')

#%%
plt.rcParams.update({'font.size': 32})
N = 100
K = 3
#clf =[LogisticRegression(), LogisticRegression()]
#clf = [RandomForestClassifier()]
#clf = [RandomForestClassifier(), LogisticRegression(), LogisticRegression(),
#       LogisticRegression(), LogisticRegression(), LogisticRegression(),
#       LogisticRegression()]
#clf= [MLPClassifier(activation = 'logistic', max_iter = 1000)]
clf = [RandomForestClassifier(), LogisticRegression(), LogisticRegression(),
       LogisticRegression(), LogisticRegression(), LogisticRegression(),
       LogisticRegression(), LogisticRegression(), LogisticRegression()]
#clf = [RandomForestClassifier(), RandomForestClassifier(), RandomForestClassifier(),
#       RandomForestClassifier(), RandomForestClassifier(), RandomForestClassifier(),
#       RandomForestClassifier(), RandomForestClassifier()]
#

#tx = [data[['ADC_med', 'max_tum_area', 'T2w_diff_var_mean', 'tum_vol',
#            'T2w_contrast_mean']],
#    data[['tum_vol', 'ADC_ave', 'T2w_ave']],
#    data[['tum_vol', 'ADC_ave', 'T2w_ave', 'total_dose']],
#    data[['8w_tum_area', '8w_tum_area_norm', 
#            '12w_tum_area', '12w_tum_area_norm', '8w_int_tum_area',
#            '8w_int_tum_area_norm', '12w_int_tum_area', 
#            '12w_int_tum_area_norm']],
#    data[['8w_tum_area_norm', '8w_int_tum_area_norm']]]

    
#tx = [data[['ADC_med', 'max_tum_area', 'T2w_diff_var_mean', 'tum_vol',
#            'T2w_contrast_mean']],
#      data[['max_tum_area', 'T2w_diff_var_mean', 'tum_vol',
#            'T2w_contrast_mean']],
#      data[['ADC_med', 'T2w_diff_var_mean', 'tum_vol', 'T2w_contrast_mean']],
#      data[['ADC_med', 'max_tum_area', 'tum_vol', 'T2w_contrast_mean']],
#      data[['ADC_med', 'max_tum_area', 'T2w_diff_var_mean',
#            'T2w_contrast_mean']],
#      data[['ADC_med', 'max_tum_area', 'T2w_diff_var_mean', 'tum_vol']]]

#
#tx = [data[['ADC_med', 'max_tum_area', 'tum_vol', 'T2w_diff_var_mean', 'T2w_contrast_mean']], 
#      data[['tum_vol', 'ADC_ave', 'T2w_ave', 'total_dose']],
#      data[['init_tum_area', 'total_dose']],
#      data[['8w_tum_area']],
#      data[['8w_int_tum_area']],
#      data[['8w_tum_area_norm']],
#      data[['8w_int_tum_area_norm']],
#      data[['12w_tum_area']],
#      data[['12w_int_tum_area']],
#      data[['12w_tum_area_norm']],
#      data[['12w_int_tum_area_norm']]]

#tx = [data[['ADC_med', 'max_tum_area', 'T2w_diff_var_mean', 'tum_vol',
#            'T2w_contrast_mean']],
#      data[['ADC_med', 'max_tum_area', 'T2w_diff_var_mean', 'tum_vol',
#            'T2w_contrast_mean', 'init_tum_area', '8w_tum_area_norm',
#            '8w_int_tum_area_norm']],
#      data[['ADC_med', 'max_tum_area', 'T2w_diff_var_mean', 'tum_vol',
#            'T2w_contrast_mean', '8w_tum_area_norm', '8w_int_tum_area_norm']],
#      data[['tum_vol', 'ADC_ave', 'T2w_ave', 'total_dose']],
#      data[['tum_vol', 'ADC_ave', 'T2w_ave', 'total_dose', 'init_tum_area',
#            '8w_tum_area', '8w_tum_area_norm']],
#      data[['tum_vol', 'ADC_ave', 'T2w_ave', 'total_dose', '8w_tum_area',
#            '8w_tum_area_norm']],
#      data[['init_tum_area', '8w_tum_area_norm', '8w_int_tum_area_norm']],
#      data[['8w_tum_area_norm', '8w_int_tum_area_norm']]]
      
tx = [data[['ADC_med', 'max_tum_area', 'T2w_diff_var_mean', 'tum_vol',
            'T2w_contrast_mean']],
      data[['1w_int_tum_area_norm']],
      data[['2w_int_tum_area_norm']],
      data[['3w_int_tum_area_norm']],
      data[['4w_int_tum_area_norm']],
      data[['5w_int_tum_area_norm']],
      data[['6w_int_tum_area_norm']],
      data[['7w_int_tum_area_norm']],
      data[['8w_int_tum_area_norm']]]


#tx = [data[['8w_tum_area', '8w_tum_area_norm', '12w_tum_area', 
#            '12w_tum_area_norm', '8w_int_tum_area', '8w_int_tum_area_norm',
#            '12w_int_tum_area', '12w_int_tum_area_norm']]]
#

#tx = [data[['init_tum_area']],
#      data[['TTum330_alphaG1120_vd03_norm']],
#      data[['TTum330_alphaG1120_norm']],
#      data[['TTum330_alphaG1120_vd05_norm']]]

#tx = [data[['ADC_med', 'max_tum_area', 'T2w_diff_var_mean', 'tum_vol',
#            'T2w_contrast_mean']],
#      data[['tum_vol', 'ADC_ave', 'T2w_ave']],
#      data[['init_tum_area']],
#      data[['TTum85']],
#      data[['TTum330']],
#      data[['mean_val']],
#      data[['TTum820']]]

#tx = [data[['8w_tum_area']],
#      data[['8w_tum_area_norm']],
#      data[['12w_tum_area']],
#      data[['12w_tum_area_norm']],
#      data[['8w_int_tum_area']],
#      data[['8w_int_tum_area_norm']],
#      data[['12w_int_tum_area']],
#      data[['12w_int_tum_area_norm']]]

y = data[['bio_rec_6']].to_numpy().ravel()  
#y = data[['gleason']].to_numpy().ravel()  

scores = np.zeros((N, K, len(tx)))
    
for j in range(N) :
        cv = StratifiedKFold(n_splits = K, shuffle = True)
        for i, x in enumerate(tx) :
            scores[j, :, i] = cross_val_score(clf[i], x, y, cv = cv,
                  scoring = 'roc_auc')
        
scoresMean = np.mean(scores, axis = 1)
#scoresMean = pd.DataFrame(data = scoresMean, columns = ['038'])

#scoresMean = pd.DataFrame(data = scoresMean, columns = ['init', '02' '03',
#                                                        '038', '05'])
      
#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['5_feat.', 'No_ADC_med',
#                                     'No_max_tum_area', 'No_T2w_diff_var_mean',
#                                     'No_tum_vol', 'No_T2w_contrast_mean'])

#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['kheamara_feat.', 
#                                     'im_feat.', 
#                                     'init_area',
#                                     '8w_tum_area',
#                                     '8w_int_tum_area',
#                                     '8w_tum_area_norm',
#                                     '8w_int_tum_area_norm',
#                                     '12w_tum_area',
#                                     '12w_int_tum_area',
#                                     '12w_tum_area_norm',
#                                     '12w_int_tum_area_norm'])  

scoresMean = pd.DataFrame(data = scoresMean,
                          columns = ['kheamara_feat.', 
                                     '1w_int_tum_area_norm', 
                                     '2w_int_tum_area_norm',
                                     '3w_int_tum_area_norm',
                                     '4w_int_tum_area_norm',
                                     '5w_int_tum_area_norm',
                                     '6w_int_tum_area_norm',
                                     '7w_int_tum_area_norm',
                                     '8w_int_tum_area_norm'])  
#          
#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['im_feat.', 'tissue_feat.',
#                                     'init_tum_area', '8w_tum_area'])        
  
#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['im_feat.', 'tissue_feat.',
#                                     'init_tum_area', '8w_tum_area'])
#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['8w_tum_area', '8w_tum_area_norm',
#                                     '12w_tum_area', '12w_tum_area_norm',
#                                     '8w_int_tum_area', '8w_int_tum_area_norm',
#                                     '12w_int_tum_area', '12w_int_tum_area_norm'])
#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['3_Khemara_feat.', 'im_feat_&_total_dose',
#                                     'init_tum_area_&_total_dose', 
#                                     '8w_tum_area_norm',
#                                     '8w_int_tum_area_norm'])
#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['260', '330', '360', '400'])
#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['khemara_feat.', 'im_feat.',
#                                     'init_tum_area', '85', '330',
#                                     '565', '820'])
#scoresMean = pd.DataFrame(data = scoresMean,
#                          columns = ['khemara_feat.', 'im_feat.',
#                                     'init_tum_area', '090', '120',
#                                     '154', '223'])

print('Mean AUC')
print(scoresMean.mean(axis = 0))

print('Standard deviation AUC')
print(scoresMean.std(axis = 0))

print('Median AUC')
print(scoresMean.median(axis = 0))
    
#%%
plt.close('all')

tcolor = [red, orange, blue, green, greenBlue, darkPurple, redOrange, redPurple,
          red, orange, blue, green, greenBlue, darkPurple, redOrange, redPurple]
ax = sns.boxplot(data = scoresMean, orient = 'v', palette = tcolor)
#ax = sns.swarmplot(data = scoresMean, color = ".25")

#add_stat_annotation(ax, data = scoresMean,
#                    box_pairs = [('init', '03'), ('03', '038')],
#                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
#                    line_offset = -0.05, verbose = 2)

#add_stat_annotation(ax, data = scoresMean,
#                    box_pairs = [('no_immuno', 'immuno')],
#                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
#                    line_offset = -0.05, verbose = 2)
#add_stat_annotation(ax, data = scoresMean,
#                    box_pairs = [('260',
#                                '330'),
#                                ('330',
#                                '400'),
#                                ('260',
#                                '400')],
#                    test = 'Wilcoxon', text_format = 'star', loc = 'outside',
#                    line_offset = -0.05, verbose = 2)
ax.set(ylabel = 'AUC', ylim = [0.6, 1])
ax.set_yticks([0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_xticklabels([], ha = 'center')

#%%   
plt.close('all')
plt.rcParams.update({'font.size': 32})  
N = 100
K = 3
#
#clf = [LogisticRegression(), LogisticRegression(), LogisticRegression(),
#       LogisticRegression()]

clf = [RandomForestClassifier(), LogisticRegression(), LogisticRegression(),
       LogisticRegression(), LogisticRegression(), LogisticRegression(),
       LogisticRegression(), LogisticRegression(), LogisticRegression()]

#tx = [data[['max_tum_area', 'ADC_ave', 'T2w_ave',]],
#      data[['noHypNec_vascDensNoPref_038_ADCT2w']],
#      data[['noHypNec_vascDensNoPref_038_simp3_ADCT2w']]]

#tx = [data[['tum_vol', 'ADC_ave', 'T2w_ave']],
#      data[['tum_area_from_vol', 'dens_ADCT2w']],
#      data[['init_tum_area']],
#      data[['TTum330_alphaG1120']]]

tx = [data[['ADC_med', 'max_tum_area', 'T2w_diff_var_mean', 'tum_vol',
            'T2w_contrast_mean']],
      data[['1w_int_tum_area_norm']],
      data[['2w_int_tum_area_norm']],
      data[['3w_int_tum_area_norm']],
      data[['4w_int_tum_area_norm']],
      data[['5w_int_tum_area_norm']],
      data[['6w_int_tum_area_norm']],
      data[['7w_int_tum_area_norm']],
      data[['8w_int_tum_area_norm']]]


y = data[['bio_rec']].to_numpy().ravel()  
#y = data[['gleason']].to_numpy().ravel()  

mean_fpr = np.linspace(0, 1, 500)
mean_tprK = np.zeros((N, 500, len(tx)))
auc_roc = np.zeros((N, len(tx))) 
auc_rec_pre = np.zeros((N, len(tx))) 

mean_recall = np.linspace(0, 1, 500)
mean_preK = np.zeros((N, 500, len(tx)))
no_skill = 0

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
        no_skillK = 0
        for k, (train, test) in enumerate(cv.split(x, y)) :
#            xS, yS = SMOTEENN().fit_sample(x[train], y[train])
            xS, yS = x[train], y[train]
            xtest, ytest = x[test], y[test]
            clf[i].fit(xS, yS)
            probas = clf[i].predict_proba(xtest)
            fpr, tpr, _ = roc_curve(ytest, probas[:, 1]) 
            mean_tprK[j, :, i] += np.interp(mean_fpr, fpr, tpr)
            mean_tprK[j, 0, i] = 0.0
            precision, recall, _ = precision_recall_curve(ytest, probas[:, 1])
            precision = np.flip(precision)
            recall = np.flip(recall)
            mean_preK[j, :, i] += np.interp(mean_recall, recall, precision)
            no_skillK += len(ytest[ytest == 1]) / len(test)
#            intercept[j, i] += clf[i].intercept_
#            coef[i][j, :] += clf[i].coef_.ravel()
            
    no_skillK /= K
    no_skill += no_skillK

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

mean_preK /= K
auc_rec_pre = np.trapz(mean_preK, mean_recall, axis = 1) 

mean_preN = np.mean(mean_preK, axis = 0)
std_preN = np.std(mean_preK, axis = 0, ddof = 1)

no_skill /= N

level = 0.95
dof = K - 1

fig, ax_roc = plt.subplots() 
#tcolor = [blue, 'tab:gray', brown]

tcolor = [red, orange, blue, green, greenBlue, darkPurple, redOrange, redPurple, redGreen]
#tcolor = [darkPurple, darkBlue, orangePurple, orange, greenBlue]

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

#ax_roc.legend(['Med. ADC, max. tum. area, tum. vol.,\nT2w diff. var. mean & '
#               'T2w contrast mean (RF)\n(AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_roc[0], std_auc_roc[0]),
#               'Tum. vol., ave. T2w, ave. ADC & d (LR)\n'
#               '(AUC = %0.2f $\pm$ %0.2f)' % (mean_auc_roc[1], std_auc_roc[1]),
#               'Init. tum. area & d (LR) (AUC = %0.2f $\pm$ %0.2f)'
#               % (mean_auc_roc[2], std_auc_roc[2]),
#               'Norm. tum. area at 8 w (LR)\n(AUC = %0.2f $\pm$ %0.2f)'
#               % (mean_auc_roc[3], std_auc_roc[3]),
#               'Norm. int. of tum. area up to 8 w (LR)'
#               '\n(AUC = %0.2f $\pm$ %0.2f)'
#               % (mean_auc_roc[4], std_auc_roc[4])], loc = 'lower right',
#             fontsize = 24)

#ax_roc.legend(['kheamara_feat. (LR) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_roc[0], std_auc_roc[0]), 
#               'kheamara_feat_init_8w_int8w (LR) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_roc[1], std_auc_roc[1]),
#               'kheamara_feat_8w_int8w (LR) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_roc[2], std_auc_roc[2]),
#               'im_feat. (LR) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_roc[3], std_auc_roc[3]), 
#               'im_feat_init_8w_int8w (LR) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_roc[4], std_auc_roc[4]),
#               'im_feat_8w_int8w (LR) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_roc[5], std_auc_roc[5]),
#               'init_8w_int8w (LR) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_roc[6], std_auc_roc[6]),
#               '8w_int8w (LR) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_roc[7], std_auc_roc[7])], loc = 'lower right',
#             fontsize = 21)

#ax_roc.legend(['kheamara_feat. (RF) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_roc[0], std_auc_roc[0]), 
#               'im_feat. (LR) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_roc[1], std_auc_roc[1]), 
#               'init. (LR) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_roc[2], std_auc_roc[2]),
#               '4w (LR) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_roc[3], std_auc_roc[3]),
#               '4w norm (LR) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_roc[4], std_auc_roc[4])], loc = 'lower right',
#             fontsize = 21)
#    
#ax_roc.legend(["Imaging parameters\n(AUC = %0.2f $\pm$ %0.2f)" %
#           (mean_auc_roc[0], std_auc_roc[0]),
#           "Tum. area at 8 w. from\ncomprehensive model\n(AUC = %0.2f $\pm$ %0.2f)" %
#           (mean_auc_roc[1], std_auc_roc[1]),
#           "Tum. area at 8 w. from\nreduced model\n(AUC = %0.2f $\pm$ %0.2f)" % 
#           (mean_auc_roc[2], std_auc_roc[2])], loc = 'lower right')

#ax_roc.legend(['Init. tum. area\n(AUC = %0.2f $\pm$ %0.2f)' %
#           (mean_auc_roc[0], std_auc_roc[0]),
#           'Norm. int. of tum. area up to 8 w.\n(vd = 0.02) (AUC = %0.2f $\pm$ %0.2f)' %
#           (mean_auc_roc[1], std_auc_roc[1]),
#           'Norm. int. of tum. area up to 8 w.\n(vd = 0.025) (AUC = %0.2f $\pm$ %0.2f)' %
#           (mean_auc_roc[2], std_auc_roc[2]),
#           'Norm. int. of tum. area up to 8 w.\n(vd = 0.03) (AUC = %0.2f $\pm$ %0.2f)' %
#           (mean_auc_roc[3], std_auc_roc[3]),
#           'Norm. int. of tum. area up to 8 w.\n(vd = 0.038) (AUC = %0.2f $\pm$ %0.2f)' % 
#           (mean_auc_roc[4], std_auc_roc[4]),
#           'Norm. int. of tum. area up to 8 w.\n(vd = 0.05) (AUC = %0.2f $\pm$ %0.2f)' % 
#           (mean_auc_roc[5], std_auc_roc[5])], loc = 'lower right')
    
#ax_roc.legend(['Khemara param.\n(AUC = %0.2f $\pm$ %0.2f)' %
#               (mean_auc_roc[0], std_auc_roc[0]),
#               'Imaging param.\n(AUC = %0.2f $\pm$ %0.2f)' %
#               (mean_auc_roc[1], std_auc_roc[1]),
#               'Init. tum. area\n(AUC = %0.2f $\pm$ %0.2f)' %
#               (mean_auc_roc[2], std_auc_roc[2]),
#               '85\n(AUC = %0.2f $\pm$ %0.2f)' %
#               (mean_auc_roc[3], std_auc_roc[3]),
#               '330\n(AUC = %0.2f $\pm$ %0.2f)' %
#               (mean_auc_roc[4], std_auc_roc[4]),
#               '565\n(AUC = %0.2f $\pm$ %0.2f)' %
#               (mean_auc_roc[5], std_auc_roc[5]),
#               '820\n(AUC = %0.2f $\pm$ %0.2f)' % 
#               (mean_auc_roc[6], std_auc_roc[6])],loc = 'lower right')

fig, ax_rec_pre = plt.subplots() 

for i in range(len(tx)) :
    ax_rec_pre.plot(mean_recall, mean_preN[:, i], linewidth = 6,
                    color = tcolor[i])
    ci = stats.t.interval(level, dof, mean_preN[:, i], std_preN[:, i] /np.sqrt(N))
    ax_rec_pre.fill_between(mean_recall, ci[0], ci[1], color = tcolor[i],
                            alpha = 0.1)
    
ax_rec_pre.plot([0, 1], [no_skill, no_skill], '--', color = 'tab:gray',
                linewidth = 6)
ax_rec_pre.set(xlabel = 'Recall', ylabel = 'Precision', ylim = [0, 1])
mean_auc_rec_pre = np.mean(auc_rec_pre, axis = 0)
std_auc_rec_pre = np.std(auc_rec_pre, axis = 0)

#ax_rec_pre.legend(['Med. ADC, max. tum. area, tum. vol.,\nT2w diff. var. mean'
#                   '& T2w contrast mean (RF)\n(AUC = %0.2f $\pm$ %0.2f)' 
#                   % (mean_auc_rec_pre[0], std_auc_rec_pre[0]),
#                   'Tum. vol., ave. T2w, ave. ADC & d (LR) '
#                   '(AUC = %0.2f $\pm$ %0.2f)'
#                   % (mean_auc_rec_pre[1], std_auc_rec_pre[1]),
#                   'Init tum. area & d (LR) (AUC = %0.2f $\pm$ %0.2f)'
#                   % (mean_auc_rec_pre[2], std_auc_rec_pre[2]),
#                   'Norm. tum. area at 8 w (LR) (AUC = %0.2f $\pm$ %0.2f)'
#                   % (mean_auc_rec_pre[3], std_auc_rec_pre[3]),
#                   'Norm. int. of tum. area up to 8 w (LR)'
#                   '(AUC = %0.2f $\pm$ %0.2f)'
#                   % (mean_auc_rec_pre[4], std_auc_rec_pre[4])],
#                 loc = 'upper right', fontsize = 21)

#ax_rec_pre.legend(['kheamara_feat. (NN) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_rec_pre[0], std_auc_rec_pre[0]), 
#               'kheamara_feat_init_8w_int8w (NN) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_rec_pre[1], std_auc_rec_pre[1]),
#               'kheamara_feat_8w_int8w (NN) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_rec_pre[2], std_auc_rec_pre[2]),
#               'im_feat. (NN) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_rec_pre[3], std_auc_rec_pre[3]), 
#               'im_feat_init_8w_int8w (NN) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_rec_pre[4], std_auc_rec_pre[4]),
#               'im_feat_8w_int8w (NN) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_rec_pre[5], std_auc_rec_pre[5]),
#               'init_8w_int8w (NN) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_rec_pre[6], std_auc_rec_pre[6]),
#               '8w_int8w (NN) (AUC = %0.2f $\pm$ %0.2f)' 
#               % (mean_auc_rec_pre[7], std_auc_rec_pre[7])], loc = 'upper right',
#             fontsize = 18)

#ax_rec_pre.legend(["Imaging parameters\n(AUC = %0.2f $\pm$ %0.2f)" %
#           (mean_auc_rec_pre[0], std_auc_rec_pre[0]),
#           "Tum. area at 8 w. from\ncomprehensive model\n(AUC = %0.2f $\pm$ %0.2f)" %
#           (mean_auc_rec_pre[1], std_auc_rec_pre[1]),
#           "Tum. area at 8 w. from\nreduced model\n(AUC = %0.2f $\pm$ %0.2f)" % 
#           (mean_auc_rec_pre[2], std_auc_rec_pre[2])], loc = 'lower right')
    
#ax.legend(['8w_tum_area (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc[0], std_auc[0]),
#           '8w_tum_area_norm (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc[1], std_auc[1]),
#           '12w_tum_area (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc[2], std_auc[2]),
#           '12w_tum_area_norm (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc[3], std_auc[3]),
#           '8w_int_tum_area (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc[4], std_auc[4]),
#           '8w_int_tum_area_norm (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc[5], std_auc[5]),
#           '12w_int_tum_area (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc[6], std_auc[6]),
#           '12w_int_tum_area_norm (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc[7], std_auc[7])],
#           loc = 'lower right')
#
#ax.legend(["Imaging parameters\n(AUC = %0.2f $\pm$ %0.2f)" %
#           (mean_auc[0], std_auc[0]),
#           "Tissue parameters\n(AUC = %0.2f $\pm$ %0.2f)" %
#           (mean_auc[1], std_auc[1]),
#           "Tum. area at 0 w.\n(AUC = %0.2f $\pm$ %0.2f)" % 
#           (mean_auc[2], std_auc[2]),
#           "Norm. int. of tum. area up to 8 w.\n(AUC = %0.2f $\pm$ %0.2f)" %
#           (mean_auc[3], std_auc[3])],
#           loc = 'lower right')

#ax.legend(["Khemara feat.\n(AUC = %0.2f $\pm$ %0.2f)" %
#           (mean_auc[0], std_auc[0]),
#           "No med. ADC\n(AUC = %0.2f $\pm$ %0.2f)" %
#           (mean_auc[1], std_auc[1]),
#           "No max. tum. area\n(AUC = %0.2f $\pm$ %0.2f)" % 
#           (mean_auc[2], std_auc[2]),
#           "No diff. var. mean T2-w\n(AUC = %0.2f $\pm$ %0.2f)" % 
#           (mean_auc[3], std_auc[3]),
#           "No tum. vol.\n(AUC = %0.2f $\pm$ %0.2f)" % 
#           (mean_auc[4], std_auc[4]),
#           "No T2-w contrast mean\n(AUC = %0.2f $\pm$ %0.2f)" % 
#           (mean_auc[5], std_auc[5])],
#           loc = 'lower right')

#ax.legend(["(AUC = %0.2f $\pm$ %0.2f)" % (mean_auc[0], std_auc[0])],
#           loc = 'lower right')

##
#_, p = stats.wilcoxon(auc_[:, 0], auc_[:, 1], alternative = 'less')
#print("4_features, comp.:", p)
#_, p = stats.wilcoxon(auc_[:, 0], auc_[:, 2], alternative = 'less')
#print("4_features, red.:", p)
#_, p = stats.wilcoxon(auc_[:, 1], auc_[:, 2])
#print("comp., red.:", p)

#%%
fig, ax_bal_acc = plt.subplots() 
tcolor = [red, orange, blue, green, greenBlue, darkPurple, redOrange,
          redPurple]
bal_acc = np.zeros((500, len(tx)))
for i in range(len(tx)) :
    bal_acc[:, i] = 0.5 * (1 + mean_tprN[:, i] - mean_fpr) 
    ax_bal_acc.plot(bal_acc[:, i], linewidth = 6, color = tcolor[i])

ax_bal_acc.set(xlabel = 'Threshold', ylabel = 'Balanced accuracy',
               ylim = [0, 1])

#ax_bal_acc.legend(['kheamara_feat. (NN)', 
#               'kheamara_feat_init_8w_int8w (NN)',
#               'kheamara_feat_8w_int8w (NN)',
#               'im_feat. (NN)', 
#               'im_feat_init_8w_int8w (NN)',
#               'im_feat_8w_int8w (NN)',
#               'init_8w_int8w (NN)',
#               '8w_int8w (NN)'],
#           loc = 'lower left', fontsize = 18)

fig, ax_f1 = plt.subplots() 
f1 = np.zeros((500, len(tx)))
for i in range(len(tx)) :
    f1[:, i] = 2.0 * (mean_recall * mean_preN[:, i]) / (mean_recall + mean_preN[:, i])
    ax_f1.plot(f1[:, i], linewidth = 6, color = tcolor[i])

ax_f1.set(xlabel = 'Threshold', ylabel = 'F1', ylim = [0, 1])

#ax_f1.legend(['kheamara_feat. (NN)', 
#               'kheamara_feat_init_8w_int8w (NN)',
#               'kheamara_feat_8w_int8w (NN)',
#               'im_feat. (NN)', 
#               'im_feat_init_8w_int8w (NN)',
#               'im_feat_8w_int8w (NN)',
#               'init_8w_int8w (NN)',
#               '8w_int8w (NN)'],
#           loc = 'upper left', fontsize = 18)

#%%
x = data['TTum330_alphaG1120'].to_numpy()

xp = np.linspace(min(x), max(x), 100)
logit = mean_intercept[3] + mean_coef[3] * xp 

fig, ax = plt.subplots()
plt.scatter(x, y)
plt.plot(xp, logit, linewidth = 6)
plt.plot(xp, 1 / (1 + np.exp(-logit)), linewidth = 6)

ax.set(title = 'LogReg', xlabel = 'TTum330_alphaG1120',
       ylabel = 'bio_rec', ylim = [-0.1, 1.1])



#%%
plt.close('all')
plt.rcParams.update({'font.size': 32})  
N = 100
K = 3

tx = [data[['ADC_med', 'max_tum_area', 'tum_vol', 'T2w_diff_var_mean',
            'T2w_contrast_mean']], 
      data[['tum_vol', 'ADC_ave', 'T2w_ave']],
      data[['init_tum_area']],
      data[['8w_tum_area_norm']],
      data[['8w_int_tum_area_norm']]]

y = data[['bio_rec_6']].to_numpy().ravel()  

n_thresholds = 100
probas = np.zeros((76, 2, len(tx)))
tn = np.zeros((N, n_thresholds, len(tx)))
fp = np.zeros((N, n_thresholds, len(tx)))
fn = np.zeros((N, n_thresholds, len(tx)))
tp = np.zeros((N, n_thresholds, len(tx)))
conf_mat = np.zeros((2, 2, n_thresholds, len(tx)))

clf = [RandomForestClassifier(), RandomForestClassifier(),
       LogisticRegression(), LogisticRegression(),  LogisticRegression()]

for j in range(N) :
#    cv = LeaveOneOut()
    cv = StratifiedKFold(n_splits = K, shuffle = True)   
    for i, x in enumerate(tx) :
        x = x.to_numpy()
        for k, (train, test) in enumerate(cv.split(x, y)) :
#            xS, yS = SMOTE().fit_sample(x[train], y[train])
            xS, yS = x[train], y[train]
            xtest, ytest = x[test], y[test]
            clf[i].fit(xS, yS)
            probas[test, :, i] += clf[i].predict_proba(xtest)
       
        t_thresholds = np.linspace(1, 0, n_thresholds)
        for l, threshold in enumerate(t_thresholds):
                ypred  = probas[:, 1, i] > threshold
                tnk, fpk, fnk, tpk = confusion_matrix(y, ypred,
                                                      labels = [0, 1]).ravel()
                tn[j, l, i] += tnk
                fp[j, l, i] += fpk
                fn[j, l, i] += fnk
                tp[j, l, i] += tpk
probas /= N

mean_tn = np.mean(tn, axis = 0)
mean_fp = np.mean(fp, axis = 0)
mean_fn = np.mean(fn, axis = 0)
mean_tp = np.mean(tp, axis = 0)

conf_mat[0, 0, :, :] = mean_tn
conf_mat[0, 1, :, :] = mean_fp
conf_mat[1, 0, :, :] = mean_fn
conf_mat[1, 1, :, :] = mean_tp

#%%
fig, ax_cal_curve = plt.subplots() 
fig, ax_cal_hist = plt.subplots() 
nbins = 5
#prob_true = np.zeros((nbins - 1, len(tx)))
#prob_pred = np.zeros((nbins - 1, len(tx)))
tcolor = [red, orange, blue, green, greenBlue, darkPurple, redOrange,
          redPurple]
for i in range(len(tx)) :
    prob_true, prob_pred =  calibration_curve(y, probas[:, 1, i].ravel(),
             n_bins = nbins, normalize = True, strategy = 'quantile')
    ax_cal_curve.plot(prob_pred, prob_true, 's-',
                      color = tcolor[i], linewidth = 6, markersize = 12)
    ax_cal_hist.hist(probas[:, 1, i], range = (0, 1), bins = nbins,
             histtype = 'step', lw = 2, color = tcolor[i], linewidth = 6)
    
ax_cal_curve.plot([0, 1], [0, 1], '--', color = 'tab:gray', linewidth = 6)

ax_cal_curve.legend(['Med. ADC, max. tum. area,\ntum. vol., '
                     'T2w diff. var. mean\n& T2w contrast mean (RF)',
                     'Tum. vol., ave. T2w &\nave. ADC (LR)',
                     'Init. tum. area(LR)',
                     'Tum. area at 8 w (LR)', 
                     'Int. of tum. area up to 8 w (LR)'],
                    loc = 'upper left', fontsize = 20)
ax_cal_curve.set(xlabel = 'Predicted risk', ylabel = 'Observed risk')

ax_cal_hist.legend(['Med. ADC, max. tum. area, tum. vol.,\n'
                     'T2w diff. var. mean & T2w contrast mean (RF)',
                     'Tum. vol., ave. T2w, ave. ADC & d (LR)',
                     'Init. tum. area & d (LR)',
                     'Tum. area at 8 w (LR)', 
                     'Int. of tum. area up to 8 w (LR)'],
                    loc = 'upper right', fontsize = 24)

ax_cal_hist.set(xlabel = 'Predicted risk', ylabel = 'Frequency')
#%%
plt.close('all')
plt.rcParams.update({'font.size': 32})  
N = 100
K = 3

tx = [data[['ADC_med', 'max_tum_area', 'tum_vol', 'T2w_diff_var_mean',
            'T2w_contrast_mean']], 
    data[['ADC_med', 'max_tum_area', 'tum_vol', 'T2w_diff_var_mean',
            'T2w_contrast_mean', 'init_tum_area']],
    data[['ADC_med', 'max_tum_area', 'tum_vol', 'T2w_diff_var_mean',
            'T2w_contrast_mean', '8w_tum_area_norm', '8w_int_tum_area_norm']],
    data[['ADC_med', 'max_tum_area', 'tum_vol', 'T2w_diff_var_mean',
            'T2w_contrast_mean', '8w_tum_area_norm']],
    data[['ADC_med', 'max_tum_area', 'tum_vol', 'T2w_diff_var_mean',
            'T2w_contrast_mean', '8w_int_tum_area_norm']],
    data[['tum_vol', 'ADC_ave', 'T2w_ave', 'total_dose']],
    data[['8w_tum_area_norm', '8w_int_tum_area_norm']]]

y = data[['gleason']].to_numpy().ravel()  

probas = np.zeros((76, 3, len(tx)))
conf_mat = np.zeros((3, 3, N, len(tx)))

#clf = [RandomForestClassifier(), RandomForestClassifier(), RandomForestClassifier(),
#       RandomForestClassifier(), RandomForestClassifier(), RandomForestClassifier(),
#       RandomForestClassifier()]

clf = [LogisticRegression(), LogisticRegression(), LogisticRegression(),
       LogisticRegression(), LogisticRegression(), LogisticRegression(),
       LogisticRegression()]

for j in range(N) :
#    cv = LeaveOneOut()
    cv = StratifiedKFold(n_splits = K, shuffle = True)   
    for i, x in enumerate(tx) :
        x = x.to_numpy()
        conf_matk = np.zeros((3, 3))
        for k, (train, test) in enumerate(cv.split(x, y)) :
            xS, yS = x[train], y[train]
            xtest, ytest = x[test], y[test]
            clf[i].fit(xS, yS)
            ypred = clf[i].predict(xtest)
#            probas[test, :, i] = clf[i].predict_proba(xtest)
            conf_matk += confusion_matrix(ytest, ypred, labels = [6, 7, 8])
            rep = classification_report(ytest, ypred, labels = [6, 7, 8])
        conf_mat[:, :, j, i] = conf_matk
mean_conf_mat = np.mean(conf_mat, axis = 2)

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