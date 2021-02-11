## Imports
import numpy as np
import cvxpy as cp

## Fitting with censored data
X_file = "~/development/cvxopt/hwk5/X.csv"
y_file = "~/development/cvxopt/hwk5/y.csv"
c_true_file = "~/development/cvxopt/hwk5/c_true.csv"

n = 20
M = 25
K = 100
D = -1.7712

X = np.loadtxt(open(X_file, "rb"), delimiter=",")
y = np.loadtxt(open(y_file, "rb"), delimiter=",")
c_true = np.loadtxt(open(c_true_file, "rb"), delimiter=",")

z = cp.Variable(K-M)
c = cp.Variable(n)
ss_low = cp.sum_squares(y - c.T @ X[:,0:M])
ss_upp = cp.sum_squares(z - c.T @ X[:,M:K])
ss_all = ss_low + ss_upp

p = cp.Problem(cp.Minimize(ss_all), [z >= D])
p.solve()

r_opt = cp.norm(c_true - c, 2) / cp.norm(c_true)
