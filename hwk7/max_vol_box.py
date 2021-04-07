import numpy as np
import cvxpy as cp

A_file = "~/development/cvxopt/hwk7/max_vol_box/A.csv"
b_file = "~/development/cvxopt/hwk7/max_vol_box/b.csv"
A = np.loadtxt(open(A_file, "rb"), delimiter=",")
b = np.loadtxt(open(b_file, "rb"), delimiter=",")

L = np.zeros_like(A)
U = np.zeros_like(A)
for i, j in it.product(range(A.shape[0]), range(A.shape[1])):
    L[i,j] = max(-A[i,j], 0)
    U[i,j] = max(A[i,j], 0)

l = cp.Variable(A.shape[1], nonneg=True)
u = cp.Variable(A.shape[1], nonneg=True)
obj = cp.geo_mean(u - l)
con = [L @ l + U @ u <= b]
p = cp.Problem(cp.Maximize(obj), con)
p.solve()

