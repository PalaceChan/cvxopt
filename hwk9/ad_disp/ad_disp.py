import numpy as np
import cvxpy as cp

n = 100; m = 30; T = 60
Acontr = np.loadtxt(open("~/development/cvxopt/hwk9/ad_disp/Acontr.csv", "rb"), delimiter=",")
Tcontr = np.loadtxt(open("~/development/cvxopt/hwk9/ad_disp/Tcontr.csv", "rb"), delimiter=",")
I = np.loadtxt(open("~/development/cvxopt/hwk9/ad_disp/I.csv", "rb"), delimiter=",")
p = np.loadtxt(open("~/development/cvxopt/hwk9/ad_disp/p.csv", "rb"), delimiter=",")
q = np.loadtxt(open("~/development/cvxopt/hwk9/ad_disp/q.csv", "rb"), delimiter=",")
R = np.loadtxt(open("~/development/cvxopt/hwk9/ad_disp/R.csv", "rb"), delimiter=",")

## Optimal N
N = cp.Variable(R.shape, nonneg = True)
s = cp.pos(q - cp.diag((Acontr.T @ N) @ Tcontr))
pen = p.T @ s
rev = cp.trace(R.T @ N)
pnl = rev - pen
cons = [np.ones(n).T @ N == I]

prb = cp.Problem(cp.Maximize(pnl), cons)
prb.solve()

## N if only show top revenue ad per period
N = np.zeros_like(R)
N[np.argmax(R, axis=0), np.arange(T)] = I
s = cp.pos(q - cp.diag((Acontr.T @ N) @ Tcontr))
pen = p.T @ s
rev = cp.trace(R.T @ N)
pnl = rev - pen
