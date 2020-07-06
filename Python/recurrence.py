import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from sklearn.metrics import auc, confusion_matrix, roc_curve, mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import KNeighborsRegressor 

plt.rcParams.update({'font.size': 32})

#%%
path = '../../Carlos/Results/40x2Gy_NoHypNec_TSim5months'
pathTissue = path + '/Tissue'
endTreat = np.loadtxt(pathTissue + '1/endTreatTumDens_0.res') 
endTreat[0] = endTreat[0] / 24
threeMon = np.loadtxt(pathTissue + '1/3MonTumDens_0.res')
threeMon[0] = threeMon[0] / 24
temp = np.loadtxt(pathTissue + '1/tumDens_0.res')

nTissues = 21
nrow, ncol = np.shape(temp)
output = np.zeros((nrow, ncol, nTissues))

for i in range(nTissues) :
    output[:, :, i] = np.loadtxt(pathTissue + str(i + 1) + '/tumDens_0.res')

output[:, 0, :] = output[:, 0, :] / 24

#%%
plt.close('all')
fig, ax = plt.subplots()
for i in range(nTissues) :
    ax.plot(output[:, 0, i], output[:, 1, i], linewidth = 6)
ax.set(xlabel = 'Time (days)', ylabel = 'Tumor density (%)',
       title = 'Tumor density')

#%%
plt.close('all')
fig, ax = plt.subplots()
for i in range(nTissues) :
    ax.plot(output[0, 1, i], output[-1, 1, i], 'o', linewidth = 6)
p = np.polyfit(output[0, 1, :], output[-1, 1, :], 1)
pts = np.linspace(min(output[0, 1, :]), max(output[0, 1, :]), 20)
ax.plot(pts, np.polyval(p, pts), '--k', linewidth = 6)
ax.set(xlabel = 'Initial tumor density', ylabel = '5 months tumor density',
       title = 'Tumor density')

t, pvalue = stats.wilcoxon(output[0, 1, :], output[-1, 1, :])

#%%
plt.close('all')
    
path = '../../Carlos/Results/Corr/Art_NoHypNec_2Gy_5Val_5Rep/'
combDens = np.loadtxt(path + 'combDens.res')
endTreat = np.loadtxt(path + 'endTreatTumDens.res') 
threeMon = np.loadtxt(path + '3MonTumDens.res')

nTissues = 100
plt.plot(threeMon[:, 0], endTreat[:, 0], 'o')

x = threeMon[:, 0].reshape(-1, 1)
y = endTreat[:, 0]
xtrain, xtest, ytrain, ytest = train_test_split(x, y)
linReg = LinearRegression()
linReg.fit(xtrain, ytrain)
ypred = linReg.predict(xtest)

fig, ax = plt.subplots()
ax.plot(ytest, ypred, 'o', linewidth = 6)
ax.set(title = 'Linear regression', xlabel = 'ytest', ylabel = 'ypred')

#%%
plt.close('all')
    
path = '../../Carlos/Results/Corr/Art_NoHypNec_2Gy_5Val_5Rep/'
combDens = np.loadtxt(path + 'combDens.res')
endTreat = np.loadtxt(path + 'endTreatTumDens.res') 
threeMon = np.loadtxt(path + '3MonTumDens.res')

nTissues = 100

x = combDens
y = endTreat[:, 0]

B = 10
nNeigh = 10
RMSE = np.zeros([nNeigh, B])

for b in range(B) :
    xtrain, xtest, ytrain, ytest = train_test_split(x, y)
    ybestPred = np.zeros(len(ytest))
    for i in range (nNeigh) :
        kNeigh = KNeighborsRegressor(n_neighbors = i + 1)
        kNeigh.fit(xtrain, ytrain)
        ypred = kNeigh.predict(xtest)
        RMSE[i, b] = np.sqrt(np.mean((ytest - ypred)**2))

fig, ax = plt.subplots()
ax.boxplot(RMSE.T)
ax.set(title = 'K-neighbors regression', xlabel = 'number of neighbors',
       ylabel = 'RMSE')


xtrain, xtest, ytrain, ytest = train_test_split(x, y)
kNeigh = KNeighborsRegressor(n_neighbors = 2)
kNeigh.fit(xtrain, ytrain)
ypred = kNeigh.predict(xtest)
 
fig, ax = plt.subplots()
ax.plot(ytest, '-o', linewidth = 4)
ax.plot(ypred, 'o', linewidth = 4)
ax.set(title = 'K-neighbors regression', xlabel = 'n', ylabel = 'Tumor density')
ax.legend(['Data', 'Prediction'])
