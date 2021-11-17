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
plt.plot(floatData[0],floatData[1])
plt.plot(floatData[0],floatData[2])
plt.plot(floatData[0],floatData[3])

doubleData = np.loadtxt("./build/2dlaplaceDouble.txt").transpose()
plt.plot(doubleData[0],doubleData[1])
plt.plot(doubleData[0],doubleData[2])
plt.plot(doubleData[0],doubleData[3])

plt.xscale("log")
plt.yscale("log")
plt.grid()
plt.show()
