import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#%%
def fsigmoid(x, a, b, c):
#    return 1.0 / (1.0 + np.exp(-a * (x - b)))
    return np.exp(-a * np.exp(-b * x + c))

#%%
#path = "../../Carlos/Results/TCP_0.8_1_0.01_1/"
path = "../../Carlos/Results/TCP_0.8_1_0.03_1/"

N = 100

tcp1Gy = np.zeros(N)
tcp2Gy = np.zeros(N)
tcp3Gy = np.zeros(N)
tcp4Gy = np.zeros(N)
tcp5Gy = np.zeros(N)
    
with open(path + "doseToControl_1Gy.res", 'r') as fTcp1Gy : 
    tcp1Gy = np.fromstring(fTcp1Gy.read(), dtype = float, sep = ' ')
with open(path + "doseToControl_2Gy.res", 'r') as fTcp2Gy : 
    tcp2Gy = np.fromstring(fTcp2Gy.read(), dtype = float, sep = ' ')
with open(path + "doseToControl_3Gy.res", 'r') as fTcp3Gy : 
    tcp3Gy = np.fromstring(fTcp3Gy.read(), dtype = float, sep = ' ')
with open(path + "doseToControl_4Gy.res", 'r') as fTcp4Gy : 
    tcp4Gy = np.fromstring(fTcp4Gy.read(), dtype = float, sep = ' ')
with open(path + "doseToControl_5Gy.res", 'r') as fTcp5Gy : 
    tcp5Gy = np.fromstring(fTcp5Gy.read(), dtype = float, sep = ' ')

ecdf1Gy = ECDF(tcp1Gy[np.nonzero(tcp1Gy)])
ecdf2Gy = ECDF(tcp2Gy[np.nonzero(tcp2Gy)])
ecdf3Gy = ECDF(tcp3Gy[np.nonzero(tcp3Gy)])
ecdf4Gy = ECDF(tcp4Gy[np.nonzero(tcp4Gy)])
ecdf5Gy = ECDF(tcp5Gy[np.nonzero(tcp5Gy)])

plt.plot(ecdf1Gy.x, ecdf1Gy.y, 'r')
plt.plot(ecdf2Gy.x, ecdf2Gy.y, 'c')
plt.plot(ecdf3Gy.x, ecdf3Gy.y, 'y')
plt.plot(ecdf4Gy.x, ecdf4Gy.y, 'b')
plt.plot(ecdf5Gy.x, ecdf5Gy.y, 'g')
plt.title("Tumour control probability")
plt.xlabel("Dose (Gy)")
plt.ylabel("TCP")
plt.grid(True)
plt.xlim(0, 100)
plt.legend(["1 Gy", "2 Gy", "3 Gy", "4 Gy", "5 Gy"])

#%%
td = [1, 2, 3, 4, 5]
tcolor = ['r', 'c', 'y', 'b', 'g']
nrow = 90
ncol = 90
tumDens = 0.8
N0 = tumDens * nrow * ncol
alpha = 0.146
beta = 0.146 / 5.5
DTotMax = 100
Ttum = 565.2
gamma = np.log(2) / Ttum
nPoints = 1000
h = DTotMax / nPoints
DTot = np.zeros(nPoints)
TCP  = np.zeros(nPoints)
for i, d in enumerate(td):
    N = DTotMax / d
    T = N // 5 * 7 * 24 + N % 5 * 24  
    for j in range(nPoints):
        DTot[j] = j * h
        TCP[j]  = np.exp(- N0 * np.exp(-alpha * DTot[j] *
           (1. + beta /alpha * d) + gamma * T))
    plt.plot(DTot, TCP, tcolor[i])
plt.legend(["1 Gy", "2 Gy", "3 Gy", "4 Gy", "5 Gy"])

#%%
plt.close('all')
#path = "../../Carlos/Results/TCP_Histo_5_NoHypNec/"
#path = "../../Carlos/Results/TCP_Histo_10/"
#path = "../../Carlos/Results/TCP_Histo_10_NoAng/"
#path = "../../Carlos/Results/TCP_Histo_10_NoRes/"
#path = "../../Carlos/Results/TCP_Histo_10_NoAngNoRes/"
path = "../../Carlos/Results/TCP_Histo_10_NoHypNec/"
#path = "../../Carlos/Results/TCP_Histo_5/"

nTissues = 6
P = 10
td = [1, 2, 3, 4, 5]
tcolor = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']

tcp = np.zeros(nTissues * P)
figTcp, axTcp = plt.subplots();

for k, eld in enumerate(td):
    for i in range(nTissues):
        pathTissue = path + "Tissue" + str(i + 1) + "/"
        with open(pathTissue + "doseToControl_" + 
              str(eld) + "Gy.res", 'r') as fTcp: 
            tcp[i * P:(i + 1) * P] = np.fromstring(fTcp.read(),
               dtype = float, sep = ' ')
    ecdf = ECDF(tcp[np.nonzero(tcp)])    
    axTcp.plot(ecdf.x, ecdf.y, color = tcolor[k])
    popt, pcov = curve_fit(fsigmoid, ecdf.x[1:], ecdf.y[1:],
                           p0 = [5000, 0.15, 3.0])
    sigmoid = np.zeros(len(ecdf.x) - 1)
    for j, x in enumerate(ecdf.x[1:]):
        sigmoid[j] = fsigmoid(x, popt[0], popt[1], popt[2])
    axTcp.plot(ecdf.x[1:], sigmoid, color = tcolor[k])
    
axTcp.set(title = "TCP", xlabel = "Total dose (Gy)", ylabel = "TCP",
          xlim = [0, 100])
axTcp.grid(True)
axTcp.legend(["1 Gy", "2 Gy", "3 Gy", "4 Gy", "5 Gy"])

#%%
td = [1, 2, 3, 4, 5]
tcolor = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
N0 = 4417
alpha = 0.146
beta = 0.146 / 5.5
DTotMax = 300
Ttum = 565.2
gamma = np.log(2) / Ttum
nPoints = 1000
h = DTotMax / nPoints
DTot = np.zeros(nPoints)
TCP  = np.zeros(nPoints)
for i, d in enumerate(td):
    N = DTotMax / d
    T = N // 5 * 7 * 24 + N % 5 * 24 
    for j in range(nPoints):
        DTot[j] = j * h
        TCP[j]  = np.exp(- N0 * np.exp(-alpha * (1. + beta / alpha * d) * 
           DTot[j] + gamma * T))
    plt.plot(DTot, TCP, tcolor[i])
plt.legend(["1 Gy", "2 Gy", "3 Gy", "4 Gy", "5 Gy"])

#%%
plt.close('all')
path = "../../Carlos/Results/"
process = ["TCP_Histo_10/", "TCP_Histo_10_NoAng/", "TCP_Histo_10_NoRes/",
           "TCP_Histo_10_NoAngNoRes/", "TCP_Histo_10_NoHypNec/"]
processNames = ["All processes", "No angiogenesis", "No resorption",
                "No angiogenesis, no resorption", "No hypoxic necrosis"]
td = [1, 2, 3, 4, 5]
nTissues = 21
P = 10
tcp = np.zeros(nTissues * P)
tcp50 = np.zeros([len(td), len(process)])
tcolor = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
tpos = [-0.2, -0.1, 0.0, 0.1, 0.2]

figTcp50, axTcp50 = plt.subplots();
    
for k, eld in enumerate(td):
    figTcp, axTcp = plt.subplots();
    for i, elp in enumerate(process):
        pathProcess = path + elp
        for j in range(nTissues):
            pathTissue = pathProcess + "Tissue" + str(j + 1) + "/"
            with open(pathTissue + "doseToControl_" +
                      str(eld) + "Gy.res", 'r') as fTcp: 
                        tcp[j * P:(j + 1) * P] = np.fromstring(fTcp.read(), 
                           dtype = float, sep = ' ')
        ecdf = ECDF(tcp[np.nonzero(tcp)])    
        axTcp.plot(ecdf.x, ecdf.y)
        tcp50[k, i] = ecdf.x[ecdf.y > 0.5][0]
        axTcp50.bar(eld + tpos[i], tcp50[k, i], width = 0.1, color = tcolor[i])
        
    axTcp.set(title = "Tumour control probability (" + str(eld) + " Gy)", 
          xlabel = "Total dose (Gy)", ylabel = "TCP", xlim = [0, 100])
    axTcp.grid(True)
    axTcp.legend(processNames)
    
    axTcp50.set(title = "TCP50", xlabel = "Dose per session (Gy)", 
               ylabel = "Total dose to TCP50 (Gy)")
    axTcp50.grid(True)
    axTcp50.set_axisbelow(True)
    axTcp50.legend(processNames)
    



# las curvas teoricas de TCP estan a la izquierda de las curvas experimentales