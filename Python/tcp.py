import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt

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
#path = "../../Carlos/Results/TCP_Histo_5_NoHypNec/"
#path = "../../Carlos/Results/TCP_Histo_10/"
#path = "../../Carlos/Results/TCP_Histo_10_NoAng/"
#path = "../../Carlos/Results/TCP_Histo_10_NoRes/"
path = "../../Carlos/Results/TCP_Histo_10_NoAngNoRes/"
#path = "../../Carlos/Results/TCP_Histo_5/"

nTissues = 21
P = 10

tcp1Gy = np.zeros(nTissues * P)
tcp2Gy = np.zeros(nTissues * P)
tcp3Gy = np.zeros(nTissues * P)
tcp4Gy = np.zeros(nTissues * P)
tcp5Gy = np.zeros(nTissues * P)

for i in range(nTissues):
    pathTissue = path + "Tissue" + str(i + 1) + "/"
    
    with open(pathTissue + "doseToControl_1Gy.res", 'r') as fTcp1Gy: 
        tcp1Gy[i * P : (i + 1) * P] = np.fromstring(fTcp1Gy.read(), 
              dtype = float, sep = ' ')
    with open(pathTissue + "doseToControl_2Gy.res", 'r') as fTcp2Gy: 
        tcp2Gy[i * P : (i + 1) * P] = np.fromstring(fTcp2Gy.read(),
              dtype = float, sep = ' ')
    with open(pathTissue + "doseToControl_3Gy.res", 'r') as fTcp3Gy: 
        tcp3Gy[i * P : (i + 1) * P] = np.fromstring(fTcp3Gy.read(),
              dtype = float, sep = ' ')
    with open(pathTissue + "doseToControl_4Gy.res", 'r') as fTcp4Gy: 
        tcp4Gy[i * P : (i + 1) * P] = np.fromstring(fTcp4Gy.read(),
              dtype = float, sep = ' ')
    with open(pathTissue + "doseToControl_5Gy.res", 'r') as fTcp5Gy: 
        tcp5Gy[i * P : (i + 1) * P] = np.fromstring(fTcp5Gy.read(),
              dtype = float, sep = ' ')

ecdf1Gy = ECDF(tcp1Gy[np.nonzero(tcp1Gy)])
ecdf2Gy = ECDF(tcp2Gy[np.nonzero(tcp2Gy)])
ecdf3Gy = ECDF(tcp3Gy[np.nonzero(tcp3Gy)])
ecdf4Gy = ECDF(tcp4Gy[np.nonzero(tcp4Gy)])
ecdf5Gy = ECDF(tcp5Gy[np.nonzero(tcp5Gy)])

plt.plot(ecdf1Gy.x, ecdf1Gy.y)
plt.plot(ecdf2Gy.x, ecdf2Gy.y)
plt.plot(ecdf3Gy.x, ecdf3Gy.y)
plt.plot(ecdf4Gy.x, ecdf4Gy.y)
plt.plot(ecdf5Gy.x, ecdf5Gy.y)
plt.title("Tumour control probability")
plt.xlabel("Dose (Gy)")
plt.ylabel("TCP")
plt.grid(True)
plt.xlim(0, 100)
plt.legend(["1 Gy", "2 Gy", "3 Gy", "4 Gy", "5 Gy"])

#%%
td = [1, 2, 3, 4, 5]
tcolor = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
N0 = 4417
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
nfig = 0;
path = "../../Carlos/Results/"
process = ["TCP_Histo_10/", "TCP_Histo_10_NoAng/", "TCP_Histo_10_NoRes/",
           "TCP_Histo_10_NoAngNoRes/", "TCP_Histo_5_NoHypNec/"]
processNames = ["All processes", "No angiogenesis", "No resorption",
                "No angiogenesis, no resorption", "No hypoxic necrosis"]
doses = [1, 2, 3, 4, 5]
nTissues = 21
P = 10
tcp = np.zeros(nTissues * P)
tcp50 = np.zeros([len(doses), len(process)])
tcolor = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
tpos = [-0.2, -0.1, 0.0, 0.1, 0.2]

figTcp50, axTcp50 = plt.subplots();
    
for k, eld in enumerate(doses):
    figTcp, axTcp = plt.subplots();
    for i, elp in enumerate(process):
        if(i == 4):
            P = 5
        else:
            P = 10
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
        
    axTcp.set(title = "Tumour control probability (" + str(k) + " Gy)", 
          xlabel = "Total dose (Gy)", ylabel = "TCP", xlim = [0, 100])
    axTcp.grid(True)
    axTcp.legend(processNames)
    
    axTcp50.set(title = "TCP50", xlabel = "Dose per session (Gy)", 
               ylabel = "Total dose to TCP50")
    axTcp50.grid(True)
    axTCP50.set_axisbelow(True)
    axTcp50.legend(processNames)
    




