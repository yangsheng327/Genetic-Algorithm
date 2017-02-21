# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 17:15:06 2017

@author: Sunrise
"""

#Genetic Algorithm
import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib import ticker, cm

class GA:
#==============================================================================
#initialize the class   
#==============================================================================
    def __init__(self, fitfun, lower, upper, bit_num, chrom_num = 32):
        self.fitfun = fitfun
        self.df = lower.shape[0]
        chrom_num = chrom_num + (chrom_num%2)
        self.chrom_num = chrom_num
        self.lower = lower
        self.upper = upper
        
        self.bit_start = np.zeros(self.df)#starting index of each parameter
        self.bit_end = np.zeros(self.df)
        for p in range(1,self.df):
            self.bit_start[p] = sum(bit_num[:p])     
        self.bit_end[:-1] = self.bit_start[1:]
        self.bit_end[-1] = sum(bit_num)#
        
        self.chroms = np.random.randint(0, 2, [chrom_num, sum(bit_num)])
        self.paras = np.zeros([chrom_num, self.df])
        self.fitness = np.zeros(chrom_num)
        self.totalfit = 0
#==============================================================================
#parse chroms to parameters and update the fitness of each chrom     
#==============================================================================
    def chroms2paras(self):
        
        def gene2scale(gene):
            x = 0
            for b in range(len(gene)):
                x = 2.*x + gene[b]
            return x/(2**len(gene)-1)
        
        for n in range(self.chrom_num):
            for p in range(self.df):
                self.paras[n,p] = self.lower[p] + \
                    (self.upper[p]-self.lower[p])* \
                        gene2scale(self.chroms[n,self.bit_start[p]:
                            self.bit_end[p]])
            self.fitness[n] = self.fitfun(self.paras[n,:])
        self.totalfit = sum(self.fitness)
#==============================================================================
#reproduction with roulette-wheel
#==============================================================================
    def reproduce(self):
        
        def rouletteWheel(cmf):
            seed = np.random.uniform(0,1)
#            print(seed)
            for i in range(len(cmf)):
                if seed<= cmf[i]:
                    return i
                                
        chroms_new = np.zeros([self.chrom_num, self.bit_end[-1]])
        cmf = np.zeros(self.chrom_num)
        for n in range(self.chrom_num):
            cmf[n] = sum(self.fitness[:n+1])/self.totalfit
#        print(cmf)
        for n in range(self.chrom_num):
#            print(rouletteWheel(cmf))
            chroms_new[n:,] = self.chroms[rouletteWheel(cmf)]
        self.chroms = chroms_new
#==============================================================================
#crossover
#==============================================================================
    def crossover(self):
        
        def switch(chrom1, chrom2, pos):
            return np.append(chrom1[:pos+1], chrom2[pos+1:]), \
                np.append(chrom2[:pos+1], chrom1[pos+1:])
        
        pair_id = np.arange(self.chrom_num)
        np.random.shuffle(pair_id)
#        print(pair_id)
        for n in range(int(self.chrom_num/2)): #pair 2n and 2n+1
            pos = np.random.randint(0,self.bit_end[-1]) #switch include [pos]
#            print('pos',pos)
            self.chroms[pair_id[2*n]], self.chroms[pair_id[2*n+1]] = \
                switch(self.chroms[pair_id[2*n]], self.chroms[pair_id[2*n+1]], pos)
        
#==============================================================================
#Mutation
#==============================================================================
    def mutate(self, mut_num):
        for i in range(mut_num):
            mut_chrom = np.random.randint(0,self.chrom_num)
            mut_pos = np.random.randint(0,self.bit_end[-1])
#            print ('chrom,pos:', mut_chrom,mut_pos)
            self.chroms[mut_chrom, mut_pos] = \
                (self.chroms[mut_chrom, mut_pos] +1)%2
                
#==============================================================================
#interation
#==============================================================================
    def iterate(self, iter_num):
        self.iter_totalfit = np.zeros(iter_num+1)
        self.chroms2paras()
        self.iter_totalfit[0] = self.totalfit
        for i in range(iter_num):
            self.reproduce()
            self.crossover()
            self.mutate(1)
            self.chroms2paras()
            self.iter_totalfit[i+1] = self.totalfit
        
#def loss(x):
#    r = np.sqrt(sum(x**2)) + np.random.uniform(0,0.5)
#    return -(np.sin(r)/r)
#            
#def fitfun(x):
#    return np.exp(-loss(x))
#
##Optimization    
#l = np.array([-20, -20])
#u = np.array([20, 20])
#bit_num = [24]*len(l)
#iter_num = 200
#    
##plot
#N = 100
#x = np.linspace(l[0], u[0], N)
#y = np.linspace(l[1], u[1], N)
#X, Y = np.meshgrid(x, y)
#Z = np.zeros([N,N])
#for i in range(N):
#    for j in range(N):
#        Z[i,j] = loss(np.append(X[i,j],Y[i,j]))
#levels = np.linspace(np.amin(Z),np.amax(Z), 100)
##plt.figure(figsize=(10,10))
#plt.contourf(X, Y, Z, levels, cmap=cm.rainbow)
#plt.colorbar()
#    
#
#
#ga = GA(fitfun, l, u, bit_num, chrom_num=128)
#ga.chroms2paras()
#plt.plot(ga.paras[:,0], ga.paras[:,1], 'ko')
#
#ga.iterate(iter_num)
#plt.plot(ga.paras[:,0], ga.paras[:,1], 'ro')
#
#
#plt.figure()
#plt.plot(ga.iter_totalfit,'-')

   
