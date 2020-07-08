# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pandas as pd
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
import seaborn as sns
from scipy.optimize import curve_fit
import scipy.stats as stats
from sklearn.metrics import auc, confusion_matrix, roc_curve
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from scipy.interpolate import LinearNDInterpolator
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import KFold
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
def fsigmoid(x, a, b):
    return 1.0 / (1.0 + np.exp(-a *(x - b)))
#    return np.exp(-a * np.exp(-b * x + c))
#%%

#data = pd.read_csv('Tumor_76_ADC.csv')
data = pd.read_csv('../../Carlos/VPIPRO/Tumor_76_ADC_T2_simp.csv')
output = np.loadtxt('../../Carlos/Results/Khemara/output.res')
data['fin_tum_area'] = output;
#data = data.loc[data['hormono'] == 0]
plt.close('all')
nfig = 0

#%%
nfig += 1
plt.figure(nfig)
data.boxplot(by = 'gleason', column = ['ADC_average'], grid = False,
             fontsize = 32)

nfig += 1
plt.figure(nfig)
data.boxplot(by = 'gleason', column = ['tumor_volume'], grid = False,
             fontsize = 32)

#%%
nfig += 1
plt.figure(nfig)
data.boxplot(by = 'bio_rec', column = ['max_tumor_area'], grid = False,
             fontsize = 32)


#%%
plt.close('all')

sigmoid = np.zeros(301)
x = list(range(301))
fig, ax = plt.subplots();

adcRec = data.loc[data['bio_rec'] == 1][['ADC_average']].to_numpy()
adcRec = adcRec.ravel()
ecdfRec = ECDF(adcRec)
popt, pcov = curve_fit(fsigmoid, ecdfRec.x[1:], ecdfRec.y[1:], p0 = [1, 50])
for j in range(301):
    sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
ax.plot(ecdfRec.x, ecdfRec.y, linewidth = 6)
#ax.plot(x, sigmoid, linewidth = 6)

adcNoRec = data.loc[data['bio_rec'] == 0][['ADC_average']].to_numpy()
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

thresholds = np.linspace(float(data[['ADC_average']].min()),
                         float(data[['ADC_average']].max()), 20)

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
tx = ['ADC_average', 'max_tumor_area', 'tumor_volume', 'T2w_average',
      'epi_dens_ADC_ave', 'epi_dens_T2w_ave', 'epi_dens_mean']
txNames = ['ADC', 'maxTumArea', 'tumVol', 'T2w', 'dens(ADC)', 'dens(T2w)',
           'dens(ADC&T2w)']

#tx = ['ADC_average', 'max_tumor_area', 'T2w_average']
#txNames = ['average ADC', 'max_tumor_area', 'average T2w']

tx = ['ADC_average', 'ADC_average_norm', 'epi_dens_norm_ADC_ave',
      'epi_dens_ADC_ave', 'T2w_average', 'T2w_average_norm',
      'epi_dens_norm_T2w_ave', 'epi_dens_T2w_ave', 'epi_dens_mean']
txNames = ['ADC', 'ADCNorm', 'densNorm(ADC)', 'dens(ADC)', 'T2w', 'T2wNorm', 
           'densNorm(T2w)', 'dens(T2w)', 'dens(ADC&T2w)']
N = 100
k = 3
scores = np.zeros((N, k, len(tx)))
logReg = LogisticRegression()

for i, el in enumerate(tx) :
#    fig, ax = plt.subplots()
##    ax.hist([data.loc[data['bio_rec'] == 1][el].to_numpy(),
##             data.loc[data['bio_rec'] == 0][el].to_numpy()])
#    sns.distplot(data.loc[data['bio_rec'] == 1][el].to_numpy())
#    sns.distplot(data.loc[data['bio_rec'] == 0][el].to_numpy())
#    ax.set(title = 'Density - ' + txNames[i], xlabel = 'txNames[i]',
#           ylabel = 'density')
#    ax.legend(['Biochemical recurrence', 'No biochemical recurrence'])
    
    x = data[[el]].to_numpy()
    y = data[['bio_rec']].to_numpy().ravel() 
    
    for j in range(N) :
        cv = StratifiedKFold(n_splits = k, shuffle = True)
        scores[j, :, i] = cross_val_score(logReg, x, y, cv = cv,
              scoring = 'roc_auc')
        
scoresMean = np.mean(scores, axis = 1)
fig, ax = plt.subplots()
ax.boxplot(scoresMean)
ax.set(title = 'AUC', ylim = [0,1], ylabel = 'AUC')
ax.set_xticklabels(txNames, rotation = 45, ha = 'right')

#%%
plt.close('all')
#tx = ['ADC_average', 'max_tumor_area', 'tumor_volume', 'T2w_average',
#      'epi_dens_ADC_ave', 'epi_dens_T2w_ave', 'epi_dens_mean']
#txNames = ['ADC', 'maxTumArea', 'tumVol', 'T2w', 'dens(ADC)', 'dens(T2w)',
#           'dens(ADC&T2w)']

#tx = ['ADC_average', 'ADC_average_norm', 'epi_dens_norm_ADC_ave',
#      'epi_dens_ADC_ave', 'T2w_average', 'T2w_average_norm',
#      'epi_dens_norm_T2w_ave', 'epi_dens_T2w_ave', 'epi_dens_mean',
#      'fin_tum_area']
#txNames = ['ADC', 'ADCNorm', 'densNorm(ADC)', 'dens(ADC)', 'T2w', 'T2wNorm', 
#           'densNorm(T2w)', 'dens(T2w)', 'dens(ADC&T2w)', 'finTumArea']

tx = ['ADC_average', 'epi_dens_ADC_ave', 'max_tumor_area', 'fin_tum_area']
txNames = ['ADC', 'dens(ADC)', 'maxTumArea','finTumArea']

#tx = ['epi_dens_ADC_ave', 'epi_dens_T2w_ave', 'epi_dens_mean', 'random']
#txNames = ['dens(ADC)', 'dens(T2w)', 'dens(ADC&T2w)', 'random']

N = 1000
K = 9
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
            fpr, tpr, _= roc_curve(y[test], probas[:, 0])   
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
    ax.set(title = 'ROC - ' + txNames[i], xlabel = 'FPR', ylabel = 'TPR',
           ylim = [0, 1])
    ax.legend(['ROC (AUC = %0.2f)' % auc_], loc = 'lower right')
    
    
x = data[['ADC_average', 'max_tumor_area']].to_numpy()
    
fig, ax = plt.subplots() 
mean_fpr = np.linspace(0, 1, 100)
mean_tprK = np.zeros((N, 100))   
    
for j in range(N) :
    cv = StratifiedKFold(n_splits = K, shuffle = True)
                 
#        fig, ax = plt.subplots() 
    for k, (train, test) in enumerate(cv.split(x, y)) :
        logReg.fit(x[train], y[train])
        probas = logReg.predict_proba(x[test])
        fpr, tpr, _= roc_curve(y[test], probas[:, 0]) 
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
ax.set(title = 'ROC - ADC and maxTumArea', xlabel = 'FPR', ylabel = 'TPR',
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