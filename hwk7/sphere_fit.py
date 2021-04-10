import numpy as np
import cvxpy as cp

U_file = "~/development/cvxopt/hwk7/sphere_fit_data/U.csv"
U = np.loadtxt(open(U_file, "rb"), delimiter=",")

x = cp.Variable(2); r = cp.Variable(1, nonneg=True)

## Solve Minimum radius sphere containing the points
obj = r
cons = [cp.norm(U[:,i] - x) / r <= 1 for i in range(50)]

p = cp.Problem(cp.Minimize(obj), cons)
p.solve(qcp=True, solver=cp.SCS)

## Knowing the center, solve for the radius
r2 = cp.Variable(1, nonneg=True)
obj = cp.norm(cp.hstack([cp.norm(U[:,i] - x.value)**2 - r2 for i in range(50)]))

p = cp.Problem(cp.Minimize(obj))
p.solve()
print(np.sqrt(r2.value))
