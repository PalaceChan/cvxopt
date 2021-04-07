import numpy as np
import cvxpy as cp

M = 20;
N = 20;
P = 20;

X_file = "~/development/cvxopt/hwk7/sep3way_data/X.csv"
Y_file = "~/development/cvxopt/hwk7/sep3way_data/Y.csv"
Z_file = "~/development/cvxopt/hwk7/sep3way_data/Z.csv"

X = np.loadtxt(open(X_file, "rb"), delimiter=",")
Y = np.loadtxt(open(Y_file, "rb"), delimiter=",")
Z = np.loadtxt(open(Z_file, "rb"), delimiter=",")

a1 = cp.Variable(2); a2 = cp.Variable(2); a3 = cp.Variable(2)
b1 = cp.Variable(1); b2 = cp.Variable(1); b3 = cp.Variable(1)
u = cp.Variable(20, nonneg=True); v = cp.Variable(20, nonneg=True); w = cp.Variable(20, nonneg=True)

obj = cp.sum(u) + cp.sum(v) + cp.sum(w)
c1s = [a1.T @ X[:,i] - b1 >= cp.max(cp.hstack([a2.T @ X[:,i] - b2, a3.T @ X[:,i] - b3])) - u[i] for i in range(20)]
c2s = [a2.T @ Y[:,i] - b2 >= cp.max(cp.hstack([a1.T @ Y[:,i] - b1, a3.T @ Y[:,i] - b3])) - v[i] for i in range(20)]
c3s = [a3.T @ Z[:,i] - b3 >= cp.max(cp.hstack([a1.T @ Z[:,i] - b1, a2.T @ Z[:,i] - b2])) - w[i] for i in range(20)]

p = cp.Problem(cp.Minimize(obj), c1s + c2s + c3s)
p.solve()
