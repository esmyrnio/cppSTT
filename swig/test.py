import pySTT
import matplotlib.pyplot as plt
import numpy as np
import math

size = 100 #number of models to be computed
minDensity = 0.5 #(x10^15cm/gr^3)
maxDensity = 4.0 #(x10^15cm/gr^3)
densities = np.logspace(math.log10(minDensity),math.log10(maxDensity),size)
theoryName = "DEF" #thoryName = "R2"
eosName = "ppwff1.cold.rns1.1.txt" #see EoS/ for other equations of state available
coupling = -5.0 #coupling = 20.0

mass,radius,scalarCharge,minimizationError,centralScalar = [np.zeros(size) for i in range(5)]
for centralDensity,i in zip(densities,range(size)):
    mass[i],radius[i],scalarCharge[i],minimizationError[i],centralScalar[i] = pySTT.get_MR(theoryName,eosName,centralDensity,coupling)
    print(centralDensity,mass[i],radius[i],scalarCharge[i],minimizationError[i],centralScalar[i])
 
# plt.plot(radius,mass)
# plt.show()
