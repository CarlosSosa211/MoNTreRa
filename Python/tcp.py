import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt

#path = "../../Carlos/Results/TCP_0.8_1_0.1_1/"
path = "../../Carlos/Results/TCP_Histo_5/"

nTissues = 1
P = 5

tcp1Gy = np.zeros(nTissues * P)
tcp2Gy = np.zeros(nTissues * P)
tcp3Gy = np.zeros(nTissues * P)
tcp4Gy = np.zeros(nTissues * P)
tcp5Gy = np.zeros(nTissues * P)

for i in range(nTissues):
    pathTissue = path + "Tissue" + str(i + 1) + "/"
    
    with open(pathTissue + "doseToControl_1Gy.res", 'r') as fTcp1Gy : 
        tcp1Gy[i * P : (i + 1) * P] = np.fromstring(fTcp1Gy.read(), dtype = float, sep = ' ')
    with open(pathTissue + "doseToControl_2Gy.res", 'r') as fTcp2Gy : 
        tcp2Gy[i * P : (i + 1) * P] = np.fromstring(fTcp2Gy.read(), dtype = float, sep = ' ')
    with open(pathTissue + "doseToControl_3Gy.res", 'r') as fTcp3Gy : 
        tcp3Gy[i * P : (i + 1) * P] = np.fromstring(fTcp3Gy.read(), dtype = float, sep = ' ')
    with open(pathTissue + "doseToControl_4Gy.res", 'r') as fTcp4Gy : 
        tcp4Gy[i * P : (i + 1) * P] = np.fromstring(fTcp4Gy.read(), dtype = float, sep = ' ')
    with open(pathTissue + "doseToControl_5Gy.res", 'r') as fTcp5Gy : 
        tcp5Gy[i * P : (i + 1) * P] = np.fromstring(fTcp5Gy.read(), dtype = float, sep = ' ')

ecdf1Gy = ECDF(tcp1Gy)
ecdf2Gy = ECDF(tcp2Gy)
ecdf3Gy = ECDF(tcp3Gy)
ecdf4Gy = ECDF(tcp4Gy)
ecdf5Gy = ECDF(tcp5Gy)

plt.plot(ecdf1Gy.x, ecdf1Gy.y, 'r')
plt.plot(ecdf2Gy.x, ecdf2Gy.y, 'c')
plt.plot(ecdf3Gy.x, ecdf3Gy.y, 'y')
plt.plot(ecdf4Gy.x, ecdf4Gy.y, 'b')
plt.plot(ecdf5Gy.x, ecdf5Gy.y, 'g')
plt.title("Tumour control probability")
plt.xlabel("Dose (Gy)")
plt.ylabel("TCP")
plt.grid(True)
plt.legend(["1 Gy", "2 Gy", "3 Gy", "4 Gy", "5 Gy"])
plt.show
