import numpy as np
import cvxpy as cp

## Relaxed LP
A = np.loadtxt(open(A_file, "rb"), delimiter=",")
b = np.loadtxt(open(b_file, "rb"), delimiter=",")
c = np.loadtxt(open(c_file, "rb"), delimiter=",")

x = cp.Variable(c.shape[0])
p = cp.Problem(cp.Minimize(c.T@x), [A @ x <= b, 0 <= x, x <= 1])
p.solve()

for t in np.linspace(0.0, 1.0, num=101):
    g = np.where(x.value > t, 1, 0)
    v = c.T @ g
    m = np.max(A @ g - b)
    if m <= 0:
        print(f"with t={t}, m={m} and feasible objective val is {v}")
        break

## Simple portfolio optimization
pbar = np.loadtxt(open(pbar_file, "rb"), delimiter=",")
S = np.loadtxt(open(S_file, "rb"), delimiter=",")
xu = np.loadtxt(open(xu_file, "rb"), delimiter=",")
wans = np.ones(xu.shape[0])

xu_ret = pbar.T @ xu
xu_rsk = np.sqrt(xu.T @ (S @ xu))

#min risk portfolio with same return as uniform
x = cp.Variable(xu.shape[0])
p = cp.Problem(cp.Minimize(cp.quad_form(x, S)), [pbar.T @ x == xu_ret, wans.T @ x == 1])
mr_rsk = np.sqrt(p.solve())

#min risk portfolio with same return as uniform but long-only
p = cp.Problem(cp.Minimize(cp.quad_form(x, S)), [pbar.T @ x == xu_ret, wans.T @ x == 1, 0 <= x])
lo_rsk = np.sqrt(p.solve())

#now limit the total short position to be <= 0.5
p = cp.Problem(cp.Minimize(cp.quad_form(x, S)), [pbar.T @ x == xu_ret, wans.T @ x == 1, wans.T @ cp.neg(x) <= 0.5])
sc_rsk = np.sqrt(p.solve())
