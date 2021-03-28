import numpy as np
import cvxpy as cp

N = 100
xtrue = np.zeros((N,1));
xtrue[0:39,] = 0.1
xtrue[49,] = 2
xtrue[69:79,] = 0.15
xtrue[79,] = 1
xtrue = xtrue.cumsum()

h = np.array([1, -0.85, 0.7, -0.3])
k = len(h)
yhat = cp.conv(h,xtrue).value

noise = np.array([-0.43,-1.7,0.13,0.29,-1.1,1.2,1.2,-0.038,0.33,0.17,-0.19,0.73,-0.59,2.2,-0.14,0.11,1.1,0.059,-0.096,-0.83,0.29,-1.3,0.71,1.6,-0.69,0.86,1.3,-1.6,-1.4,0.57,-0.4,0.69,0.82,0.71,1.3,0.67,1.2,-1.2,-0.02,-0.16,-1.6,0.26,-1.1,1.4,-0.81,0.53,0.22,-0.92,-2.2,-0.059,-1,0.61,0.51,1.7,0.59,-0.64,0.38,-1,-0.02,-0.048,4.3e-05,-0.32,1.1,-1.9,0.43,0.9,0.73,0.58,0.04,0.68,0.57,-0.26,-0.38,-0.3,-1.5,-0.23,0.12,0.31,1.4,-0.35,0.62,0.8,0.94,-0.99,0.21,0.24,-1,-0.74,1.1,-0.13,0.39,0.088,-0.64,-0.56,0.44,-0.95,0.78,0.57,-0.82,-0.27])
y = yhat[0:N] + noise

A = np.zeros(shape=(N,N))
for i in range(n):
    for j in range(len(h)):
        if i - j >= 0:
            A[i, i - j] = h[j]

x = cp.Variable(N)
fn = cp.norm(A @ x - y)
cNN = x >= 0
cMO = cp.diff(x) >= 0

p_ml = cp.Problem(cp.Minimize(fn), [cNN, cMO])
p_ml.solve()
x_ml = x.value

p_mlf = cp.Problem(cp.Minimize(fn))
p_mlf.solve()
x_mlf = x.value

#Plotting via R
#",".join([str(k) for k in xtrue]) x3 then
#x <- rbindlist(list(data.table(x=1:100, y = xtrue, type = "tru"), data.table(x=1:100, y=xml, type = "ml"), data.table(x=1:100, y = xmlf, type = "mlf")))
#ggplot(x, aes(y=y, x=x, col=type)) + geom_point() + geom_line()
