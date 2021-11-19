import numpy as np
import matplotlib.pyplot as plt

data32 = np.loadtxt("./build/laplace32.txt").transpose()
plt.plot(data32[0], data32[1], label="Single precision, order 4")
plt.plot(data32[0], data32[2], label="Single precision, order 6")
data64 = np.loadtxt("./build/laplace64.txt").transpose()
plt.plot(data64[0], data64[1], label="Double precision, order 4")
plt.plot(data64[0], data64[2], label="Double precision, order 6")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("$N$")
plt.ylabel("error")
plt.grid()
plt.legend()
plt.savefig("./document/figs/1d.pdf")
plt.show()