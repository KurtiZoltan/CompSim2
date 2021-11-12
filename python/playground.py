import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg.decomp_svd import null_space

def factorial(i):
    if i == 0:
        return 1
    else:
        return i * factorial(i-1)

def taylorPolinom(hs, order):
    derivatives = [0] * len(hs)
    ret = [[1, derivatives]]
    for direction, h in enumerate(hs):
        temp = []
        for i in range(order+1):
            derivatives = [0] * len(hs)
            derivatives[direction] = i
            temp.append([h**i / factorial(i), derivatives])
        newRet = []
        for a in ret:
            for b in temp:
                newRet.append([a[0] * b[0], [x + y for x, y in zip(a[1], b[1])]])
        ret = []
        for a in newRet:
            found = False
            for b in ret:
                if a[1] == b[1]:
                    b[0] += a[0]
                    found = True
            if not found and sum(a[1]) <= order:
                ret.append(a)
    return ret

def derivatives(variableCount, order):
    derivatives = [0] * variableCount
    ret = [derivatives]
    for direction in range(variableCount):
        temp = []
        for i in range(order+1):
            derivatives = [0] * variableCount
            derivatives[direction] = i
            temp.append(derivatives)
        newRet = []
        for a in ret:
            for b in temp:
                newRet.append([x + y for x, y in zip(a, b)])
        ret = []
        for a in newRet:
            found = False
            for b in ret:
                if a == b:
                    found = True
            if not found and sum(a) <= order:
                ret.append(a)
    return ret

def kernelPoints(radius, dim):
    if dim == 1:
        return [[i] for i in range(radius + 1)]
    else:
        ret = []
        ending = kernelPoints(radius, dim - 1)
        for i in range(radius + 1):
            for end in ending:
                if i <= end[0]:
                    new = end.copy()
                    new.insert(0, i)
                    ret.append(new)
    return ret

def extendBySymmetry(positions, symmetryGenerator):
    for p in positions:
        if symmetryGenerator(p) not in positions:
            positions.append(symmetryGenerator(p))

def mirror(p, coordinate):
    p[coordinate] *= -1
    return p

def swap(p, i, j):
    temp = p[i]
    p[i] = p[j]
    p[j] = temp
    return p

def kernelNDMatrix(radius, order, dim):
    columns = []
    for position in kernelPoints(radius, dim):
        positionList = [position]
        for i in range(dim):
            extendBySymmetry(positionList, lambda p: mirror(p.copy(), i))
        for i in range(dim):
            for j in range(i, dim):
                extendBySymmetry(positionList, lambda p: swap(p.copy(), i, j))
                
        taylors = [taylorPolinom(p, order) for p in positionList]
        taylorSum = [0] * len(taylors[0])
        for t in taylors:
            taylorSum = [x + y[0] for x, y in zip(taylorSum, t)]
        columns.append(taylorSum)
    return np.array(columns).transpose()

def kernel2DMatrix(radius, order):
    return kernelNDMatrix(radius, order, 2)

def kernel3DMatrix(radius, order):
    return kernelNDMatrix(radius, order, 3)

def kernelFromCoeffs (a, radius, dim):
    ret = np.zeros(shape=[2*radius+1]*dim)
    points = kernelPoints(radius, dim)
    for idx in np.ndindex(ret.shape):
        key = np.array(idx)
        key = key - radius
        key = np.abs(key)
        key = sorted(key)
        key = list(key)
        ret[idx] = a[points.index(key)]
    return ret

def laplace(order, dim):
    basis = derivatives(dim, order)
    ret = np.zeros((len(basis)))
    for i, b in enumerate(basis):
        for j in range(dim):
            secondDiff = [0] * dim
            secondDiff[j] = 2
            if b == secondDiff:
                ret[i] = 1
    return ret

def laplaceKernelsND(radius, dim):
    order = 2 * radius
    A = kernelNDMatrix(radius, order, dim)
    b = laplace(order, dim)
    x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
    nullSpace = null_space(A)
    return x, nullSpace



def smallestError2D(radius):
    order = radius * 2
    A = kernel2DMatrix(radius, order)
    basis = derivatives(2, order)
    b = np.zeros((len(basis)))
    for i, diff in enumerate(basis):
        if diff == [2, 0] or diff == [0, 2]:
            b[i] = 1

    nullspace = null_space(A)
    a0, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)

    #print(f"Nullspace:\n", nullspace)
    #print(f"Particular solution:\n", a0)

    result = A @ a0
    #print("Result: ", list(filter(lambda x: abs(x[0]) > 1e-6, [[x, y] for x, y in zip(result, basis)])))

    #a0 + nullspace @ x are the solutions
    A = kernel2DMatrix(radius, order + 4)
    basis = derivatives(2, order + 4)
    b = np.zeros((len(basis)))
    for i, diff in enumerate(basis):
        if diff == [2, 0] or diff == [0, 2]:
            b[i] = 1

    #A @ (a0 + nullspace @ x) =~ b
    #A @ nullspace @ x =~ b - A @ a0
    x, residuals, rank, s = np.linalg.lstsq(A @ nullspace, b - A @ a0, rcond=None)
    print(f"Final result for parameters:\n", a0 + nullspace @ x)
    result = A @ (a0 + nullspace @ x)
    print(f"Final result for estimate:\n", list(filter(lambda x: abs(x[0]) > 1e-6, [[x, y] for x, y in zip(result, basis)])))
    return a0 + nullspace @ x

def drawKesrnel2D(a, radius):
    coeffs = np.zeros((2*radius+1, 2*radius+1))
    for i in range(-radius, radius + 1):
        for j in range(-radius, radius + 1):
            i0 = abs(i)
            j0 = abs(j)
            if i0 > j0:
                temp = i0
                i0 = j0
                j0 = temp
            coeffs[i + radius, j + radius] = a[i0 * (radius + 1 + radius + 1 - (i0 - 1)) // 2 + j0 - i0]
    plt.imshow(coeffs)
    plt.colorbar()
    plt.show()

radius = 2
#a = smallestError2D(radius)
#drawKernel2D(a, radius)
x, nullspace = laplaceKernelsND(radius, 10)
print(nullspace.shape)
