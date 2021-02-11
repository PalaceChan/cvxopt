## Imports
import numpy as np
import cvxpy as cp

k = 201
t = np.array([-3 + 6*(i-1)/(k-1) for i in range(1,k+1)])
y = np.exp(t)

a = cp.Variable(3)
b = cp.Variable(2)
z = cp.Parameter()

phi1 = cp.maximum(*[(a[0] + a[1]*t[i] + a[2]*t[i]**2) - (1 + b[0]*t[i] + b[1]*t[i]**2)*(y[i] + z) for i in range(k)]) <= 0
phi2 = cp.maximum(*[(1 + b[0]*t[i] + b[1]*t[i]**2)*(y[i] - z) - (a[0] + a[1]*t[i] + a[2]*t[i]**2) for i in range(k)]) <= 0
p = cp.Problem(cp.Minimize(1), [phi1, phi2])

l = 0
u = np.exp(3)
e = 1e-3
for _ in range(100):
    z.value = (l + u) / 2
    print(f"[l, u] = [{l}, {u}], err={u-l}")
    try:
        s = p.solve()
    except cp.SolverError:
        s = p.solve(solver=cp.SCS)
    if p.status == 'infeasible':
        l = z.value
    else:
        u = z.value
        if u - l <= e:
            break
