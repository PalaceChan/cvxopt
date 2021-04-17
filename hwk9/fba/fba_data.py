import numpy as np
import cvxpy as cp
import itertools as it

## Data and problem setup
S = np.array([[1,0,0,0,0,0],[-1,1,0,0,0,0],[-1,0,1,0,0,0],[0,-1,0,2,-1,0],[0,0,0,0,1,0],[0,-2,1,0,0,1],[0,0,-1,1,0,0],[0,0,0,0,0,-1],[0,0,0,-1,0,0]]).T
orig_max = np.array([10.10,100,5.90,100,3.70,100,100,100,100])
vmax = cp.Parameter(9)
v = cp.Variable(vmax.shape[0], nonneg = True)
cons = [S @ v == 0, v <= vmax]

## Solve first part
vmax.value = orig_max
p = cp.Problem(cp.Maximize(v[8]), cons)
p.solve()
gmin = 0.2*p.solve() #save 0.2*optimal value as gmin for later

## Essentials and lethals
def is_essential(i):
    vmax.value = orig_max.copy()
    vmax.value[i] = 0
    grate = p.solve()
    return grate < gmin

def is_lethal(i, j):
    vmax.value = orig_max.copy()
    vmax.value[i] = 0
    vmax.value[j] = 0
    grate = p.solve()
    return grate < gmin    

for i in range(vmax.shape[0]):
    if is_essential(i):
        print(f"gene {i+1} was essential as grate is now {grate}")

for i, j in it.product(range(vmax.shape[0]), repeat=2):
    if i < j:
        if not is_essential(i) and not is_essential(j) and is_lethal(i, j):
            print(f"gene pair G{i+1}, G{j+1} is lethal")
            
