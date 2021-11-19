import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg.decomp_svd import null_space

def multiply(p1, p2, order):
    ret = []
    for a in p1:
        for b in p2:
            newTerm = [a[0] * b[0], [x + y for x, y in zip(a[1], b[1])]]
            found = False
            for r in ret:
                if r[1] == newTerm[1]:
                    found = True
                    r[0] += newTerm[0]
            if not found and sum(newTerm[1]) <= order:
                ret.append(newTerm)
    return ret

def factorial(i):
    if i == 0:
        return 1
    else:
        return i * factorial(i-1)

def taylorPolinom(hs, order):
    '''
    returns the taylor polinum up to order at the position displaced by hs
    '''
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
    '''
    returns the deriaties up to order listed in the same way as returned by taylorPolinom
    '''
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
    '''
    up to mirroring returns the independent points of the stencil
    '''
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
    '''
    appends to position the symmetry transformed elements of points
    '''
    for p in positions:
        if symmetryGenerator(p) not in positions:
            positions.append(symmetryGenerator(p))

def mirror(p, coordinate):
    '''
    mirrors the coordinate coordinate of point p
    '''
    p[coordinate] *= -1
    return p

def swap(p, i, j):
    '''
    swaps the ith and jth component of p
    '''
    temp = p[i]
    p[i] = p[j]
    p[j] = temp
    return p

def kernelNDMatrix(radius, order, dim):
    '''
    returns matrix A, where A @ coefficients of the stencil = coefficients of the taylor expansion
    '''
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
    '''
    returns the full stencil constructed from the independent elements of the stencil
    '''
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
    """
    returns the derivatives corresponding to the laplace operator in the basis used by taylorPolinom
    """
    basis = derivatives(dim, order)
    ret = np.zeros((len(basis)))
    for i, b in enumerate(basis):
        for j in range(dim):
            secondDiff = [0] * dim
            secondDiff[j] = 2
            if b == secondDiff:
                ret[i] = 1
    return ret

def laplaceKernelND(radius, dim):
    '''
    returns a particular solution for the independent coefficients of the stencil, and the nullspace
    nullspace @ x does not contribute derivaties up to order 2 * radius + 1, it represents the freedom left to choose the stencil
    '''
    order = 2 * radius
    A = kernelNDMatrix(radius, order, dim)
    b = laplace(order, dim)
    x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=1e-6)
    nullSpace = null_space(A, rcond=1e-6)
    return x, nullSpace

def symmetricErrorLaplaceKernelND(radius, dim):
    '''
    same as laplaceKernelND with the additional constraint that the leading order error term will be proportional to a power of the laplace operator, meaning the error will be spherically symmetric in leading order.
    '''
    order = 2 * radius + 2
    x, nullspace = laplaceKernelND(radius, dim)
    A = kernelNDMatrix(radius, order, dim)
    basis = derivatives(dim, order)
    lap = laplace(order, dim)
    a = [[x, y] for x, y in zip(lap, basis)]
    error = [[0, basis] for basis in basis]
    error[0][0] = 1
    for i in range(radius + 1):
        error = multiply(error, a, order)
    error = np.array([x[0] for x in error])
    P = np.eye(len(error)) - np.outer(error, error) / np.dot(error, error)
    #A @ (x + nullspace @ t) - b = alpha * b ^ (order + 2)
    #A @ x - b + A @ nullspace @ t = alpha * b ^ (order + 2)
    #P @ A @ x - P @ b + P @ A @ nullspace @ t = 0
    sphericalNullspace = null_space(P @ A @ nullspace, rcond=1e-6)
    t0, residuals, rank, s = np.linalg.lstsq(P @ A @ nullspace, P @ lap - P @ A @ x, rcond=1e-6)
    finalx = x + nullspace @ t0
    finalNullspace = nullspace @ sphericalNullspace
    nonsphericalError = P @ (A @ (x + nullspace @ t0) - lap)
    nonsphericalError = [[x, y] for x, y in zip(nonsphericalError, basis)]
    print(f"Non spherical error:")
    for term in nonsphericalError:
        if abs(term[0]) > 1e-10:
            print(term)
    return finalx, finalNullspace

def smallestError2D(radius):
    '''
    not used in the final report
    '''
    order = radius * 2
    A = kernel2DMatrix(radius, order)
    basis = derivatives(2, order)
    b = np.zeros((len(basis)))
    for i, diff in enumerate(basis):
        if diff == [2, 0] or diff == [0, 2]:
            b[i] = 1

    nullspace = null_space(A, rcond=1e-6)
    a0, residuals, rank, s = np.linalg.lstsq(A, b, rcond=1e-6)

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

radius = 2
dim = 2
x, null = symmetricErrorLaplaceKernelND(radius, dim)
null = null.flatten()
x += -x[-1] / null[-1] * null # only use this if radius = 2 dim = 2 and you want to obtain a stencil that has 0 in the corners
print(f"x:\n", x * 60)
print(f"Null:\n", null)
print(kernelFromCoeffs(x, radius, dim))
plt.imshow(kernelFromCoeffs(x, radius, dim))
plt.colorbar()
plt.show()
