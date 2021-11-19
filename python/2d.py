import numpy as np
import matplotlib.pyplot as plt

'''exact = np.loadtxt("./build/exact.txt").transpose()
numerica = np.loadtxt("./build/numerica.txt").transpose()
numericb = np.loadtxt("./build/numericb.txt").transpose()
numericc = np.loadtxt("./build/numericc.txt").transpose()

diff = np.abs(numerica - exact)[2:-3, 2:-3]
print("a:", np.sum(diff) / np.prod(diff.shape))
plt.imshow(np.abs(numerica - exact)[2:-3, 2:-3], origin="upper")
plt.colorbar()
plt.show()

diff = np.abs(numericb - exact)[2:-3, 2:-3]
print("b:", np.sum(diff) / np.prod(diff.shape))
plt.imshow(np.abs(numericb - exact)[2:-3, 2:-3], origin="upper")
plt.colorbar()
plt.show()

diff = np.abs(numericc - exact)[2:-3, 2:-3]
print("c:", np.sum(diff) / np.prod(diff.shape))
plt.imshow(np.abs(numericc - exact)[2:-3, 2:-3], origin="upper")
plt.colorbar()
plt.show()'''

floatData = np.loadtxt("./build/2dlaplaceFloat.txt").transpose()
plt.plot(floatData[0],floatData[1], label="Single precision, isotropic 3x3")
plt.plot(floatData[0],floatData[2], label="Single precision, anisotropic 3x3")
plt.plot(floatData[0],floatData[3], label="Single precision, isotropic 5x5")

doubleData = np.loadtxt("./build/2dlaplaceDouble.txt").transpose()
plt.plot(doubleData[0],doubleData[1], label="Double precision, isotropic 3x3")
plt.plot(doubleData[0],doubleData[2], label="Double precision, anisotropic 3x3")
plt.plot(doubleData[0],doubleData[3], label="Double precision, isotropic 5x5")

plt.xscale("log")
plt.yscale("log")
plt.xlabel("$N$")
plt.ylabel("error")
plt.grid()
plt.legend()
plt.savefig("./document/figs/2d.pdf")
plt.show()
