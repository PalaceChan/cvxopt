import numpy as np
import cvxpy as cp
import itertools as it

u1, u2 = 8.0, 20.0
s1, s2 = 6.0, 17.5
rho = -0.25

n = 100
R = np.linspace(-30, 70, num=n)

den1 = sum([np.exp(-(r - u1)**2 / (2*s1**2)) for r in R])
den2 = sum([np.exp(-(r - u2)**2 / (2*s2**2)) for r in R])
p1 = np.array([np.exp(-(r - u1)**2 / (2*s1**2)) / den1 for r in R])
p2 = np.array([np.exp(-(r - u2)**2 / (2*s2**2)) / den2 for r in R])

Y = np.zeros(shape=(n,n))
Z = np.zeros(shape=(n,n))
for i, j in it.product(range(n), repeat=2):
    Z[i,j] = (R[i] - u1)*(R[j] - u2)
    if R[i] + R[j] <= 0:
        Y[i,j] = 1


X = cp.Variable(shape=(n,n))
fn = cp.trace(X.T @ Y)
cR1 = (np.ones(n).T @ X == p2)
cR2 = (X @ np.ones(n) == p1)
cCOV = (cp.trace(Z.T @ X) == rho*s1*s2)
cPos = (X >= 0)

p = cp.Problem(cp.Maximize(fn), [cR1, cR2, cCOV, cPos])
p.solve(solver=cp.SCS)
