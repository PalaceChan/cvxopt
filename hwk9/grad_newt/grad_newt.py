import numpy as np
import cvxpy as cp

A_file = "~/development/cvxopt/hwk9/A.csv"
A = np.loadtxt(open(A_file, "rb"), delimiter=",")

x = cp.Variable(A.shape[1])
obj0 = -cp.sum([cp.log(1 - A[i,:] @ x) for i in range(A.shape[0])])
obj1 = -cp.sum([cp.log(1 - x[i]**2) for i in range(A.shape[1])])
obj = obj0 + obj1
cons = [A @ x <= 1, cp.abs(x) <= 1]

p = cp.Problem(cp.Minimize(obj), cons)
p.solve()
