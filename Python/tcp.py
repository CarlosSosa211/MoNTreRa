import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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
lightPurple = [146/255, 141/255, 251/255]
redOrange = [255/255, 149/255, 114/255]
redPurple = [207/255, 122/255, 162/255]
orangePurple = [207/255, 147/255, 148/255]
brown = [168/255, 120/255, 110/255]
#%%
plt.rcParams.update({'font.size': 32})
#path = "../../Carlos/Results/TCP/0.8_1_0.01_1/"
path = "../../Carlos/Results/TCP/0.8_1_0.03_1/"

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
plt.title("Tumor control probability")
plt.xlabel("Dose (Gy)")
plt.ylabel("TCP")
plt.grid(True)
plt.xlim(0, 100)
plt.legend(["1 Gy", "2 Gy", "3 Gy", "4 Gy", "5 Gy"])

#%%
plt.rcParams.update({'font.size': 22})
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
plt.rcParams.update({'font.size': 32})
path = "../../Carlos/Results/TCP/Histo_10/"
#path = "../../Carlos/Results/TCP/Histo_10_NoAng/"
#path = "../../Carlos/Results/TCP/Histo_10_NoRes/"
#path = "../../Carlos/Results/TCP/Histo_10_NoArrest/"
#path = "../../Carlos/Results/TCP/Histo_10_NoHypNec/"
#path = "../../Carlos/Results/TCP/Histo_10_NoOxy/"
#path = "../../Carlos/Results/TCP/Histo_10_NoAngNoRes/"
#path = "../../Carlos/Results/TCP/Histo_10_NoAngNoArrest/"
#path = "../../Carlos/Results/TCP/Histo_10_NoResNoArrest/"
#path = "../../Carlos/Results/TCP/Histo_10_NoAngNoResNoArrest/"
#path = "../../Carlos/Results/TCP/Histo_10_NoAngNoHypNec/"


nTissues = 21
P = 10
#td = [1, 2, 3, 4, 5]
td = [2]
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
    #axTcp.plot(ecdf.x, ecdf.y, color = tcolor[k])
    popt, pcov = curve_fit(fsigmoid, ecdf.x[1:], ecdf.y[1:], p0 = [1, 50])
    print(eld, "Gy", popt)
    #sigmoid = np.zeros(len(ecdf.x) - 1)
    #for j, x in enumerate(ecdf.x[1:]):
     #   sigmoid[j] = fsigmoid(x, popt[0], popt[1])
    #axTcp.plot(ecdf.x[1:], sigmoid, color = tcolor[k])
    
    sigmoid = np.zeros(101)
    x = list(range(101))
    for j in range(101):
        sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
    axTcp.plot(x, sigmoid, color = tcolor[k], linewidth = 6)
    axTcp.set(title = "TCP (2 Gy)", xlabel = "Total dose (Gy)", ylabel = "TCP",
          xlim = [0, 100])
axTcp.grid(True)
#axTcp.legend(["1 Gy", "2 Gy", "3 Gy", "4 Gy", "5 Gy"])

#%%
plt.rcParams.update({'font.size': 32})
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
        TCP[j]  = np.exp(- N0 * np.exp(-alpha * (1. + beta / alpha * d) * 
           DTot[j] + gamma * T))
    plt.plot(DTot, TCP, tcolor[i])
plt.legend(["1 Gy", "2 Gy", "3 Gy", "4 Gy", "5 Gy"])

#%%
plt.close('all')
plt.rcParams.update({'font.size': 32})
path = "../../Carlos/Results/"
process = ["TCP/Histo_10/", "TCP/Histo_10_NoAng/", "TCP/Histo_10_NoRes/",
           "TCP/Histo_10_NoArrest/", "TCP/Histo_10_NoHypNec/"]
processNames = ["All processes", "No angiogenesis", "No healthy cell division",
                "No cycle arrest", "No hypoxic necrosis"]
td = [1, 2, 3, 4, 5]
nTissues = 21
P = 10
tcp = np.zeros(nTissues * P)
tcp50 = np.zeros([len(td), len(process)])
tcolor = ['tab:gray', red, orange, darkPurple, darkGreen]
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
        #axTcp.plot(ecdf.x, ecdf.y, color = tcolor[i], linewidth = 2)
        popt, pcov = curve_fit(fsigmoid, ecdf.x[1:], ecdf.y[1:], p0 = [1, 50])
        sigmoid = np.zeros(101)
        x = list(range(101))
        for j in range(101):
            sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
        axTcp.plot(x, sigmoid, color = tcolor[i], linewidth = 6)
        tcp50[k, i] = ecdf.x[ecdf.y > 0.5][0]
        axTcp50.bar(eld + tpos[i], tcp50[k, i], width = 0.1, color = tcolor[i])
        
    axTcp.set(title = "Tumor control probability (" + str(eld) + " Gy)", 
          xlabel = "Total dose (Gy)", ylabel = "TCP", xlim = [0, 100])
    axTcp.grid(True)
    axTcp.legend(processNames)
    figTcp.savefig("TCP1Pro" + str(eld) + "Gy.pdf", bbox_inches="tight")
    
    axTcp50.set(title = "TCP50", xlabel = "Dose per session (Gy)", 
               ylabel = "Total dose to TCP50 (Gy)")
    axTcp50.grid(True)
    axTcp50.set_axisbelow(True)
    axTcp50.legend(processNames)
    
#%%
plt.close('all')
plt.rcParams.update({'font.size': 32})
path = "../../Carlos/Results/"
process = ["TCP/Histo_10/", "TCP/Histo_10_NoAngNoRes/", 
           "TCP/Histo_10_NoAngNoArrest/", "TCP/Histo_10_NoResNoArrest/"]
processNames = ["All processes", "No angiogenesis,\nno healthy cell division",
                "No angiogenesis,\nno cycle arrest",
                "No healthy cell division,\nno cycle arrest"]
td = [1, 2, 3, 4, 5]
nTissues = 21
P = 10
tcp = np.zeros(nTissues * P)
tcp50 = np.zeros([len(td), len(process)])
tcolor = ['tab:gray', redOrange, redPurple, orangePurple]
tpos = [-0.15, -0.05, 0.05, 0.15]

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
        #axTcp.plot(ecdf.x, ecdf.y, color = tcolor[i], linewidth = 2)
        popt, pcov = curve_fit(fsigmoid, ecdf.x[1:], ecdf.y[1:], p0 = [1, 50])
        sigmoid = np.zeros(101)
        x = list(range(101))
        for j in range(101):
            sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
        axTcp.plot(x, sigmoid, color = tcolor[i], linewidth = 6)
        tcp50[k, i] = ecdf.x[ecdf.y > 0.5][0]
        axTcp50.bar(eld + tpos[i], tcp50[k, i], width = 0.1, color = tcolor[i])
        
    axTcp.set(title = "Tumor control probability (" + str(eld) + " Gy)", 
          xlabel = "Total dose (Gy)", ylabel = "TCP", xlim = [0, 100])
    axTcp.grid(True)
    axTcp.legend(processNames, loc = 'upper left')
    
    axTcp50.set(title = "TCP50", xlabel = "Dose per session (Gy)", 
               ylabel = "Total dose to TCP50 (Gy)")
    axTcp50.grid(True)
    axTcp50.set_axisbelow(True)
    axTcp50.legend(processNames)
    
#%%
plt.close('all')
plt.rcParams.update({'font.size': 32})
path = "../../Carlos/Results/"
process = ["TCP/Histo_10/", "TCP/Histo_10_NoAngNoResNoArrest/"]
processNames = ["All processes",
                "No angiogenesis,\nno healthy cell division,\nno cycle arrest"]
td = [1, 2, 3, 4, 5]
nTissues = 21
P = 10
tcp = np.zeros(nTissues * P)
tcp50 = np.zeros([len(td), len(process)])
tcolor = ['tab:gray', brown]
tpos = [-0.05, 0.05]

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
        #axTcp.plot(ecdf.x, ecdf.y, color = tcolor[i], linewidth = 2)
        popt, pcov = curve_fit(fsigmoid, ecdf.x[1:], ecdf.y[1:], p0 = [1, 50])
        sigmoid = np.zeros(101)
        x = list(range(101))
        for j in range(101):
            sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
        axTcp.plot(x, sigmoid, color = tcolor[i], linewidth = 6)
        tcp50[k, i] = ecdf.x[ecdf.y > 0.5][0]
        axTcp50.bar(eld + tpos[i], tcp50[k, i], width = 0.1, color = tcolor[i])
        
    axTcp.set(title = "Tumor control probability (" + str(eld) + " Gy)", 
          xlabel = "Total dose (Gy)", ylabel = "TCP", xlim = [0, 100])
    axTcp.grid(True)
    axTcp.legend(processNames, loc = 'upper left')
    
    axTcp50.set(title = "TCP50", xlabel = "Dose per session (Gy)", 
               ylabel = "Total dose to TCP50 (Gy)")
    axTcp50.grid(True)
    axTcp50.set_axisbelow(True)
    axTcp50.legend(processNames)
        
#%%
plt.close('all')
plt.rcParams.update({'font.size': 22})
path = "../../Carlos/Results/"
process = ["TCP/Histo_10/", "TCP/Histo_10_NoAngNoHypNec/"]
processNames = ["All processes", "No angiogenesis, no hypoxic necrosis"]
td = [1, 2, 3, 4, 5]
nTissues = 21
P = 10
tcp = np.zeros(nTissues * P)
tcp50 = np.zeros([len(td), len(process)])
tcolor = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray']
tpos = [0.0, 0.1, 0.2, 0]

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
        
    axTcp.set(title = "Tumor control probability (" + str(eld) + " Gy)", 
          xlabel = "Total dose (Gy)", ylabel = "TCP", xlim = [0, 100])
    axTcp.grid(True)
    axTcp.legend(processNames)
    
    axTcp50.set(title = "TCP50", xlabel = "Dose per session (Gy)", 
               ylabel = "Total dose to TCP50 (Gy)")
    axTcp50.grid(True)
    axTcp50.set_axisbelow(True)
    axTcp50.legend(processNames)
    
#%%
plt.close('all')
plt.rcParams.update({'font.size': 32})
path = "../../Carlos/Results/"
process = ["TCP/Histo_10/", "TCP/Histo_10_NoAng/", "TCP/Histo_10_NoRes/",
           "TCP/Histo_10_NoAngNoRes/", "TCP/Histo_10_NoAngNoResNoArrest/",
           "TCP/Histo_10_NoHypNec/", "TCP/Histo_10_NoAngNoHypNec/",
           "TCP/Histo_10_NoOxy/"]
processNames = ["All processes", "No angiogenesis", "No healthy cell division",
                "No angiogenesis, no healthy cell division",
                "No angiogenesis, no resoprtion, no cycle arrest",
                "No hypoxic necrosis", "No angiogenesis, no hypoxic necrosis",
                "No oxygenation (no angiogenesis, no hypoxic necrosis)"]
td = [1, 2, 3, 4, 5]
nTissues = 21
P = 10
tcp = np.zeros(nTissues * P)
tcp50 = np.zeros([len(td), len(process)])
tcolor = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray']
tpos = [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4]

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
        
    axTcp.set(title = "Tumor control probability (" + str(eld) + " Gy)", 
          xlabel = "Total dose (Gy)", ylabel = "TCP", xlim = [0, 100])
    axTcp.grid(True)
    axTcp.legend(processNames)
    
    axTcp50.set(title = "TCP50", xlabel = "Dose per session (Gy)", 
               ylabel = "Total dose to TCP50 (Gy)")
    axTcp50.grid(True)
    axTcp50.set_axisbelow(True)
    axTcp50.legend(processNames)

#%%
plt.close('all')
plt.rcParams.update({'font.size': 22})
path = "../../Carlos/Results/"
process = ["TCP/Histo_10/", "TCP/Histo_10_TimeConstantOxy/", 
           "TCP/Histo_10_NoOxy/", "TCP/Histo_10_meanpO2/"]
processNames = ["All processes", "Time constant oxygenation",
                "No oxygenation (no angiogenesis)",
                "Mean pO2"]
#td = [1, 2, 3, 4, 5]
td = [2]
nTissues = 21
P = 10
tcp = np.zeros(nTissues * P)
tcp50 = np.zeros([len(td), len(process)])
tcolor = ['tab:blue', 'tab:red', 'tab:orange', 'tab:purple', 'tab:green']
tpos = [-0.1, 0.0, 0.1, 0.2]

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
#        axTcp.plot(ecdf.x, ecdf.y, color = tcolor[i])
        popt, pcov = curve_fit(fsigmoid, ecdf.x[1:], ecdf.y[1:], p0 = [1, 50])
        print(popt)
        sigmoid = np.zeros(101)
        x = list(range(101))
        for j in range(101):
            sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
        axTcp.plot(x, sigmoid, color = tcolor[i], linewidth = 4)
        tcp50[k, i] = ecdf.x[ecdf.y > 0.5][0]
        axTcp50.bar(eld + tpos[i], tcp50[k, i], width = 0.1, color = tcolor[i])
        
    axTcp.set(title = "Tumor control probability (" + str(eld) + " Gy)", 
    xlabel = "Total dose (Gy)", ylabel = "TCP", xlim = [0, 100])
    axTcp.grid(True)
    axTcp.legend(processNames)
    
    axTcp50.set(title = "TCP50", xlabel = "Dose per session (Gy)", 
               ylabel = "Total dose to TCP50 (Gy)")
    axTcp50.grid(True)
    axTcp50.set_axisbelow(True)
    axTcp50.legend(processNames)    

#%%
plt.close('all')
plt.rcParams.update({'font.size': 32})
path = "../../Carlos/Results/"
process = ["TCP/Histo_10_TTum_85/", "TCP/Histo_10_TTum_330/",
           "TCP/Histo_10_TTum_575/","TCP/Histo_10_TTum_820/",
           "TCP/Histo_10_TTum_1065/","TCP/Histo_10_TTum_1310/"]
processNames = ["$T_{tum}$ = 85 h", "$T_{tum}$ = 330 h", "$T_{tum}$ = 575 h",
                "$T_{tum}$ = 820 h", "$T_{tum}$ = 1065 h",
                "$T_{tum}$ = 1310 h"]
td = [2]
nTissues = 21
P = 10
tcp = np.zeros(nTissues * P)
tcp50 = np.zeros([len(td), len(process)])
tcolor = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray']
tpos = [-0.25, -0.15, -0.05, 0.05, 0.15, 0.25]

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
        #axTcp.plot(ecdf.x, ecdf.y, color = tcolor[i], linewidth = 2)
        popt, pcov = curve_fit(fsigmoid, ecdf.x[1:], ecdf.y[1:], p0 = [1, 50])
        print(popt)
        sigmoid = np.zeros(101)
        x = list(range(101))
        for j in range(101):
            sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
        axTcp.plot(x, sigmoid, color = tcolor[i], linewidth = 6)
        tcp50[k, i] = ecdf.x[ecdf.y > 0.5][0]
        axTcp50.bar(eld + tpos[i], tcp50[k, i], width = 0.1, color = tcolor[i])
        
    axTcp.set(title = "Tumor control probability (" + str(eld) + " Gy)", 
          xlabel = "Total dose (Gy)", ylabel = "TCP", xlim = [0, 100])
    axTcp.grid(True)
    axTcp.legend(processNames)
    
    axTcp50.set(title = "TCP50", xlabel = "Dose per session (Gy)", 
               ylabel = "Total dose to TCP50 (Gy)")
    axTcp50.grid(True)
    axTcp50.set_axisbelow(True)
    axTcp50.legend(processNames)
    
#%%
plt.close('all')
plt.rcParams.update({'font.size': 32})
path = "../../Carlos/Results/"
process = ["TCP/Histo_10_alphaG1_0.024/", "TCP/Histo_10_alphaG1_0.090/",
           "TCP/Histo_10_alphaG1_0.157/","TCP/Histo_10_alphaG1_0.223/",
           "TCP/Histo_10_alphaG1_0.290/","TCP/Histo_10_alphaG1_0.356/"]
processNames = ["$\\alpha_{tumG1}$ = 0.024 Gy$^{-1}$",
                "$\\alpha_{tumG1}$ = 0.090 Gy$^{-1}$",
                "$\\alpha_{tumG1}$ = 0.157 Gy$^{-1}$",
                "$\\alpha_{tumG1}$ = 0.223 Gy$^{-1}$",
                "$\\alpha_{tumG1}$ = 0.290 Gy$^{-1}$", 
                "$\\alpha_{tumG1}$ = 0.356 Gy$^{-1}$"]
td = [2]
nTissues = 21
P = 10
tcp = np.zeros(nTissues * P)
tcp50 = np.zeros([len(td), len(process)])
tcolor = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray']
tpos = [-0.25, -0.15, -0.05, 0.05, 0.15, 0.25]

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
#        axTcp.plot(ecdf.x, ecdf.y, color = tcolor[i], linewidth = 2)
        popt, pcov = curve_fit(fsigmoid, ecdf.x[1:], ecdf.y[1:], p0 = [1, 75])
        print(popt)
        sigmoid = np.zeros(101)
        x = list(range(101))
        for j in range(101):
            sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
        axTcp.plot(x, sigmoid, color = tcolor[i], linewidth = 6)
        tcp50[k, i] = ecdf.x[ecdf.y > 0.5][0]
        axTcp50.bar(eld + tpos[i], tcp50[k, i], width = 0.1, color = tcolor[i])
        
    axTcp.set(title = "Tumor control probability (" + str(eld) + " Gy)", 
          xlabel = "Total dose (Gy)", ylabel = "TCP", xlim = [0, 100])
    axTcp.grid(True)
    axTcp.legend(processNames)
    
    axTcp50.set(title = "TCP50", xlabel = "Dose per session (Gy)", 
               ylabel = "Total dose to TCP50 (Gy)")
    axTcp50.grid(True)
    axTcp50.set_axisbelow(True)
    axTcp50.legend(processNames)
    
#%%
plt.close('all')
plt.rcParams.update({'font.size': 32})
path = "../../Carlos/Results/"
process = ["TCP/Histo_10_alphaG2_0.025/", "TCP/Histo_10_alphaG2_0.096/",
           "TCP/Histo_10_alphaG2_0.167/","TCP/Histo_10_alphaG2_0.239/",
           "TCP/Histo_10_alphaG2_0.310/","TCP/Histo_10_alphaG2_0.381/"]
processNames = ["$\\alpha_{tumG2}$ = 0.025 Gy$^{-1}$",
                "$\\alpha_{tumG2}$ = 0.096 Gy$^{-1}$",
                "$\\alpha_{tumG2}$ = 0.167 Gy$^{-1}$",
                "$\\alpha_{tumG2}$ = 0.239 Gy$^{-1}$",
                "$\\alpha_{tumG2}$ = 0.310 Gy$^{-1}$", 
                "$\\alpha_{tumG2}$ = 0.381 Gy$^{-1}$"]
td = [2]
nTissues = 21
P = 10
tcp = np.zeros(nTissues * P)
tcp50 = np.zeros([len(td), len(process)])
tcolor = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray']
tpos = [-0.25, -0.15, -0.05, 0.05, 0.15, 0.25]

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
#        axTcp.plot(ecdf.x, ecdf.y, color = tcolor[i], linewidth = 2)
        popt, pcov = curve_fit(fsigmoid, ecdf.x[1:], ecdf.y[1:], p0 = [1, 60])
        print(popt)
        sigmoid = np.zeros(101)
        x = list(range(101))
        for j in range(101):
            sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
        axTcp.plot(x, sigmoid, color = tcolor[i], linewidth = 6)
        tcp50[k, i] = ecdf.x[ecdf.y > 0.5][0]
        axTcp50.bar(eld + tpos[i], tcp50[k, i], width = 0.1, color = tcolor[i])
        
    axTcp.set(title = "Tumor control probability (" + str(eld) + " Gy)", 
          xlabel = "Total dose (Gy)", ylabel = "TCP", xlim = [0, 100])
    axTcp.grid(True)
    axTcp.legend(processNames)
    
    axTcp50.set(title = "TCP50", xlabel = "Dose per session (Gy)", 
               ylabel = "Total dose to TCP50 (Gy)")
    axTcp50.grid(True)
    axTcp50.set_axisbelow(True)
    axTcp50.legend(processNames)
    
#%%
plt.close('all')
plt.rcParams.update({'font.size': 32})
path = "../../Carlos/Results/"
process = ["TCP/Histo_10_NoHypNec/","TCP/Histo_10_pO2Nec_0.26/", 
           "TCP/Histo_10_pO2Nec_0.52/","TCP/Histo_10_pO2Nec_0.78/", 
           "TCP/Histo_10_pO2Nec_1.04/","TCP/Histo_10_pO2Nec_1.3/"]
processNames = ["$pO_2^{nec}$ = 0 mmHg", "$pO_2^{nec}$ = 0.26 mmHg",
                "$pO_2^{nec}$ = 0.52 mmHg", "$pO_2^{nec}$ = 0.78 mmHg",
                "$pO_2^{nec}$ = 1.04 mmHg", "$pO_2^{nec}$ = 1.3 mmHg"]
td = [2]
nTissues = 21
P = 10
tcp = np.zeros(nTissues * P)
tcp50 = np.zeros([len(td), len(process)])
tcolor = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray']
tpos = [-0.25, -0.15, -0.05, 0.05, 0.15, 0.25]

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
#        axTcp.plot(ecdf.x, ecdf.y, color = tcolor[i], linewidth = 2)
        popt, pcov = curve_fit(fsigmoid, ecdf.x[1:], ecdf.y[1:], p0 = [1, 75])
        print(popt)
        sigmoid = np.zeros(101)
        x = list(range(101))
        for j in range(101):
            sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
        axTcp.plot(x, sigmoid, color = tcolor[i], linewidth = 6)
        tcp50[k, i] = ecdf.x[ecdf.y > 0.5][0]
        axTcp50.bar(eld + tpos[i], tcp50[k, i], width = 0.1, color = tcolor[i])
        
    axTcp.set(title = "Tumor control probability (" + str(eld) + " Gy)", 
          xlabel = "Total dose (Gy)", ylabel = "TCP", xlim = [0, 100])
    axTcp.grid(True)
    axTcp.legend(processNames)
    
    axTcp50.set(title = "TCP50", xlabel = "Dose per session (Gy)", 
               ylabel = "Total dose to TCP50 (Gy)")
    axTcp50.grid(True)
    axTcp50.set_axisbelow(True)
    axTcp50.legend(processNames)
#%%
plt.close('all')
plt.rcParams.update({'font.size': 22})
path = "../../Carlos/Results/"
process = ["TCP/Histo_10/", "TCP/Histo_10_OneAlphaBetaAllPhases/",
           "TCP/Histo_10_OneAlphaBetaAllCells/"]
processNames = ["Cell-phase-and-type dependent alpha and beta", 
                "Cell-type dependent alpha and beta",
                "Constant alpha and beta"]
td = [2]
nTissues = 21
P = 10
tcp = np.zeros(nTissues * P)
tcp50 = np.zeros([len(td), len(process)])
tcolor = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray']
tpos = [-0.1, 0.0, 0.1]

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
        #axTcp.plot(ecdf.x, ecdf.y, color = tcolor[i], linewidth = 2)
        popt, pcov = curve_fit(fsigmoid, ecdf.x[1:], ecdf.y[1:], p0 = [1, 50])
        print(popt)
        sigmoid = np.zeros(101)
        x = list(range(101))
        for j in range(101):
            sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
        axTcp.plot(x, sigmoid, color = tcolor[i], linewidth = 4)
        tcp50[k, i] = ecdf.x[ecdf.y > 0.5][0]
        axTcp50.bar(eld + tpos[i], tcp50[k, i], width = 0.1, color = tcolor[i])
        
    axTcp.set(title = "Tumor control probability (" + str(eld) + " Gy)", 
          xlabel = "Total dose (Gy)", ylabel = "TCP", xlim = [0, 100])
    axTcp.grid(True)
    axTcp.legend(processNames)
    
    axTcp50.set(title = "TCP50", xlabel = "Dose per session (Gy)", 
               ylabel = "Total dose to TCP50 (Gy)")
    axTcp50.grid(True)
    axTcp50.set_axisbelow(True)
    axTcp50.legend(processNames)
    
#%%
plt.close('all')
plt.rcParams.update({'font.size': 22})
path = "../../Carlos/Results/"
process = ["TCP/Histo_10/", "TCP/Histo_10_NoAngNoResNoArrest/",
           "TCP/Histo_10_OneAlphaBetaAllPhases/",
           "TCP/Histo_10_NoAngNoResNoArrestOneAlphaBetaAllPhases/"]
processNames = ["All processes",
                "No angiogenesis, no healthy cell division,\nno cycle arrest",
                "Same alpha and beta for all phases",
                "No angiogenesis, no healthy cell division,\nno cycle arrest, " 
                "same alpha and beta for\nall phases"]
td = [2]
nTissues = 21
P = 10
tcp = np.zeros(nTissues * P)
tcp50 = np.zeros([len(td), len(process)])
tcolor = ['tab:blue', 'tab:red', 'tab:orange', 'tab:green', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray']
tpos = [-0.15, -0.05, 0.05, 0.15]

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
        #axTcp.plot(ecdf.x, ecdf.y, color = tcolor[i], linewidth = 2)
        popt, pcov = curve_fit(fsigmoid, ecdf.x[1:], ecdf.y[1:], p0 = [1, 50])
        print(popt)
        sigmoid = np.zeros(101)
        x = list(range(101))
        for j in range(101):
            sigmoid[j] = fsigmoid(x[j], popt[0], popt[1])
        axTcp.plot(x, sigmoid, color = tcolor[i], linewidth = 4)
        tcp50[k, i] = ecdf.x[ecdf.y > 0.5][0]
        axTcp50.bar(eld + tpos[i], tcp50[k, i], width = 0.1, color = tcolor[i])
        
    axTcp.set(title = "Tumor control probability (" + str(eld) + " Gy)", 
          xlabel = "Total dose (Gy)", ylabel = "TCP", xlim = [0, 100])
    axTcp.grid(True)
    axTcp.legend(processNames)
    
    axTcp50.set(title = "TCP50", xlabel = "Dose per session (Gy)", 
               ylabel = "Total dose to TCP50 (Gy)")
    axTcp50.grid(True)
    axTcp50.set_axisbelow(True)
    axTcp50.legend(processNames)