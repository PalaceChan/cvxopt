import numpy as np
import cvxpy as cp

X_file = "~/development/cvxopt/hwk7/quad_metric_data_norng/X.csv"
Y_file = "~/development/cvxopt/hwk7/quad_metric_data_norng/Y.csv"
d_file = "~/development/cvxopt/hwk7/quad_metric_data_norng/d.csv"
X_test_file = "~/development/cvxopt/hwk7/quad_metric_data_norng/X_test.csv"
Y_test_file = "~/development/cvxopt/hwk7/quad_metric_data_norng/Y_test.csv"
d_test_file = "~/development/cvxopt/hwk7/quad_metric_data_norng/d_test.csv"

X = np.loadtxt(open(X_file, "rb"), delimiter=",")
Y = np.loadtxt(open(Y_file, "rb"), delimiter=",")
d = np.loadtxt(open(d_file, "rb"), delimiter=",")
X_test = np.loadtxt(open(X_test_file, "rb"), delimiter=",")
Y_test = np.loadtxt(open(Y_test_file, "rb"), delimiter=",")
d_test = np.loadtxt(open(d_test_file, "rb"), delimiter=",")

N = X.shape[1]
Z = X - Y
P = cp.Variable((5, 5), symmetric=True)
t1 = 1/N * (d**2).sum()
t2 = 1/N * cp.sum([ -2 * d[i] * cp.quad_form(Z[:,i], P)**(1/2) for i in range(N)])
t3 = 1/N * cp.sum([cp.quad_form(Z[:,i], P) for i in range(100)])
obj = t1 + t2 + t3

p = cp.Problem(cp.Minimize(obj), [P >> 0])
p.solve()

Z_test = X_test - Y_test
N = X_test.shape[1]
errors = [d_test[i] - (cp.quad_form(Z_test[:,i], P.value)**(1/2)).value for i in range(N)]
mse = np.sum(np.array(errors)**2) / N
