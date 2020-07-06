# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
import seaborn as sns
from scipy.optimize import curve_fit
import scipy.stats as stats
from sklearn.metrics import auc, confusion_matrix, roc_curve
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
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
data = pd.read_csv('../../Carlos/VPIPRO/Tumor_76_ADC_T2.csv')
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
ax.legend(['Biological recurrence', 'No biological recurrence'])

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
tx = ['ADC_average', 'max_tumor_area', 'tumor_volume']
txNames = ['average ADC', 'max_tumor_area', 'tumor volume']
for i, el in enumerate(tx) :
    fig, ax = plt.subplots()
#    ax.hist([data.loc[data['bio_rec'] == 1][el].to_numpy(),
#             data.loc[data['bio_rec'] == 0][el].to_numpy()])
    sns.distplot(data.loc[data['bio_rec'] == 1][el].to_numpy())
    sns.distplot(data.loc[data['bio_rec'] == 0][el].to_numpy())
    ax.set(title = 'Density - ' + txNames[i], xlabel = 'txNames[i]',
           ylabel = 'density')
    ax.legend(['Biological recurrence', 'No biological recurrence'])
    
    x = data[[el]].to_numpy()
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
    ax.set(title = 'ROC - ' + txNames[i])
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

#Hay pocos casos con recidiva