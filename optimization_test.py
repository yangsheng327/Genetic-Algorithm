# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 15:18:39 2017

@author: Sunrise
"""

#optimization test
import GA
import PSO

import numpy as np
import random
from scipy.optimize import minimize
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np

def loss(x):
    r = np.sqrt(sum(x**2)) + np.random.normal(0,1)
    return -(np.sin(r)/r)
            
def fitfun(x):
    return np.exp(-loss(x))

#Optimization    
l = np.array([-20, -20])
u = np.array([20, 20])
#Scipy minimize\
for i in range(5):
    x0 = [0, 0]
    x0[:] = [np.random.uniform(l[:], u[:])]
    res = minimize(loss, x0, method='BFGS')
    print(res.x)
#GA
bit_num = [24]*len(l)
iter_num = 200
#PSO
coe = np.array([0.1,0.7,0.5])
part_num=50   

#plot
N = 100
x = np.linspace(l[0], u[0], N)
y = np.linspace(l[1], u[1], N)
X, Y = np.meshgrid(x, y)
Z = np.zeros([N,N])
for i in range(N):
    for j in range(N):
        Z[i,j] = loss(np.append(X[i,j],Y[i,j]))
levels = np.linspace(np.amin(Z),np.amax(Z), 100)
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)




plt.figure()
plt.contourf(X, Y, Z, levels, cmap=cm.rainbow)
plt.colorbar()
pso = PSO.PSO(loss,part_num)
pso.part_init(l,u)
ga = GA.GA(fitfun, l, u, bit_num, chrom_num=128)
ga.chroms2paras()
plt.plot(ga.paras[:,0], ga.paras[:,1], 'bo')
plt.plot(pso.current_pos[:,0],pso.current_pos[:,1],'b*')

ga.iterate(iter_num)
pso.fit(coe,iter_num)
plt.plot(ga.paras[:,0], ga.paras[:,1], 'wo')
plt.plot(pso.current_pos[:,0],pso.current_pos[:,1],'w*')



