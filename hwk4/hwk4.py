## Imports
import numpy as np
import cvxpy as cp
import math, itertools as it

## Numerical perturbation analysis
# Problem setup
x = cp.Variable(2)
P = np.array([[1, -1/2], [-1/2, 2]])
q = np.array([-1, 0])
u1, u2 = cp.Parameter(), cp.Parameter()
f1 = np.array([1, 2])
f2 = np.array([1, -4])
f3 = np.array([1, 1])

constraints = [f1.T @ x <= u1, f2.T @ x <= u2, f3.T @ x >= -5]
prob = cp.Problem(cp.Minimize(cp.quad_form(x, P) + q.T @ x), constraints)

# Solve instance
u1.value, u2.value= -2, -3
val_opt = prob.solve()

x_opt = x.value
y_opt = [c.dual_value for c in constraints]

# Verify KKT conditions
kkt1 = [f1.T @ x_opt - u1.value, f2.T @ x_opt - u2.value, -f3.T @ x_opt - 5]
assert all([x <= 0 for x in kkt1]), "all primal conditions hold"
assert all([y >= 0 for y in y_opt]), "all duals nonnegative"
assert all([yf == 0 for yf in np.array(y_opt).T * kkt1]), "complementary slackness"

grad_f0_xopt = np.array([2*x_opt[0] - x_opt[1] - 1, 4*x_opt[1] - x_opt[0]])
grad_fi_xopt = np.array([[1, 2],[1, -4],[-1, -1]]).T @ np.array(y_opt)
assert (grad_f0_xopt + grad_fi_xopt == 0).all(), "gradient vanishes"

# Solve for a perturbed grid with forecasts
prob.backward()
i, min_i, min_diff = 1, None, math.inf
for d1, d2 in list(it.product([0, -.1, .1], repeat=2)):
    u1.value = -2 + d1
    u2.value = -3 + d2
    val_pred = val_opt - (np.array([d1, d2]).T @ y_opt[0:2])
    val_xact = prob.solve()
    val_diff = val_xact - val_pred
    if val_diff > 1e-5 and val_diff < min_diff:
        min_i = i
        min_diff = val_diff
    i += 1
    print(f"d1={d1}, d2={d2}, pred={val_pred}, xact={val_xact}, diff={val_diff}")
print(f"min i = {min_i}")    

## Option Price Bounds
s0, rfa, flr, cap = 1, 1.05, 0.9, 1.15
def get_collar_price(s1):
    if s1 > cap:
        return cap - s0
    elif flr <= s1 and s1 <= cap:
        return s1 - s0
    else:
        return flr - s0

vs = []
for s1 in np.linspace(0.5, 2.0, num=200):
    v = [get_collar_price(s1), s1, rfa, max(0, s1 - 1.1), max(0, s1 - 1.2), max(0, 0.8 - s1), max(0, 0.7 - s1)]
    vs.append(v)
V = np.array(vs)

x = cp.Variable()
p = np.array([1, 1, 0.06, 0.03, 0.02, 0.01])
y = cp.Variable(200)
c0 = V.T[0,:] @ y == x
c1 = (V.T[1:,:] @ y == p)
c2 = (y >= 0)
prob_min = cp.Problem(cp.Minimize(x), [c0, c1, c2])
prob_max = cp.Problem(cp.Maximize(x), [c0, c1, c2])
print(f"{prob_min.solve()} - {prob_max.solve()}")
