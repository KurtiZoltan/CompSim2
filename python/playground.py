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

def extendBySymmetry(positions, symmetryGenerator):
    for p in positions:
        if symmetryGenerator(p) not in positions:
            positions.append(symmetryGenerator(p))

def kernel2D(radius, order):
    columns = []
    for i in range(radius+1):
        for j in range(i, radius+1):
            positionList = [[i, j]]
            extendBySymmetry(positionList, lambda p: [-p[0], p[1]])
            extendBySymmetry(positionList, lambda p: [p[0], -p[1]])
            extendBySymmetry(positionList, lambda p: [p[1], p[0]])
            
            taylors = [taylorPolinom(p, order) for p in positionList]
            taylor = [0] * len(taylors[0])
            for t in taylors:
                taylor = [x + y[0] for x, y in zip(taylor, t)]
            columns.append(taylor)
    return np.array(columns).transpose()

def kernel3D(radius, order):
    columns = []
    for i in range(radius+1):
        for j in range(i, radius+1):
            for k in range(j, radius+1):
                positionList = [[i, j, k]]
                extendBySymmetry(positionList, lambda p: [-p[0], p[1], p[2]])
                extendBySymmetry(positionList, lambda p: [p[0], -p[1], p[2]])
                extendBySymmetry(positionList, lambda p: [p[0], p[1], -p[2]])
                extendBySymmetry(positionList, lambda p: [p[0], p[2], p[1]])
                extendBySymmetry(positionList, lambda p: [p[1], p[0], p[2]])
                extendBySymmetry(positionList, lambda p: [p[1], p[2], p[0]])
                extendBySymmetry(positionList, lambda p: [p[2], p[0], p[1]])
                extendBySymmetry(positionList, lambda p: [p[2], p[1], p[0]])
                
                taylors = [taylorPolinom(p, order) for p in positionList]
                taylor = [0] * len(taylors[0])
                for t in taylors:
                    taylor = [x + y[0] for x, y in zip(taylor, t)]
                columns.append(taylor)
    return np.array(columns).transpose()

def smallestError2D(radius):
    order = radius * 2
    A = kernel2D(radius, order)
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
    A = kernel2D(radius, order + 4)
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

def drawKernel2D(a, radius):
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

radius = 3
a = smallestError2D(radius)
drawKernel2D(a, radius)
