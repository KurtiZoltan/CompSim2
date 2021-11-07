import numpy as np
import matplotlib.pyplot as plt

data64 = np.loadtxt("../build/laplace64.txt").transpose()
plt.plot(data64[0], data64[1])
plt.plot(data64[0], data64[2])
data32 = np.loadtxt("../build/laplace32.txt").transpose()
plt.plot(data32[0], data32[1])
plt.plot(data32[0], data32[2])
plt.xscale("log")
plt.yscale("log")
#plt.plot(data32[0], data32[0]**(-2))
#plt.plot(data32[0], data32[0]**(-4))
plt.grid()
plt.show()