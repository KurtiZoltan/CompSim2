import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("./build/gaussianCoreT.txt").transpose()

plt.plot(data[0], data[1], label="Core temperature")
plt.xlabel("$t$")
plt.ylabel("$T_{core}$")
plt.grid()
plt.legend()
plt.savefig("./document/figs/gauss.pdf")
plt.show()

data = np.loadtxt("./build/conductorCoreT.txt").transpose()

plt.plot(data[0], data[1], label="Core temperature")
plt.xlabel("$t$")
plt.ylabel("$T_{core}$")
plt.grid()
plt.legend()
plt.savefig("./document/figs/conductor.pdf")
plt.show()