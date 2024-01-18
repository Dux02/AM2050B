# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 17:14:42 2024

@author: Blas
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
import time
starttime = time.time()



def KeplerEq(E, e, M):
    return E - e*np.sin(E) - M    

def KeplerEqDfn(E,e):
    return 1 - e*np.cos(E)

def KeplerSolver(e,M):
    accuracy = 1e-8
    maxIterations = 100
    if (e > 0.8):
        E = np.pi
    else:
        E = M
    
    i = 1
    while ( i < maxIterations):
        nextVal = E - (KeplerEq(E,e,M) / KeplerEqDfn(E,e))
        difference = abs(E - nextVal)
        E = nextVal
        if (difference < accuracy):
            break
        i += 1
    return E

class KeplerObject:
    
    def __init__(self,a,e,i,Omega,omega,L,m,m0,T):
        self.a = a
        self.e = e
        R = Rotation.from_euler("ZXZ", [-omega, -i, -Omega])
        self.rotMatrix = R
        self.T = T
        self.M = L - omega - Omega
        self.M0 = m0
        self.m = m
        E = KeplerSolver(e,self.M)
        nu = 2 * np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))
        r = a*(1-e*np.cos(E))
        self.pos = [r*np.cos(nu), r*np.sin(nu),0]
        self.pos = self.rotMatrix.apply(self.pos)
        
    def update(self,dt):
        self.M += dt * 2 * np.pi / self.T
        E = KeplerSolver(self.e,self.M)
        #print(E)
        nu = 2 * np.arctan(np.sqrt((1+self.e)/(1-self.e))*np.tan(E/2))
        r = self.a*(1-self.e*np.cos(E))
        self.pos = [r*np.cos(nu), r*np.sin(nu),0]
        self.pos = self.rotMatrix.apply(self.pos)
        

N = 200000
data = [[] for i in range(N)]
DT = 3600
AU = 149597870700
m0 = 1.989e30
DTS = 24*60*60
Earth = KeplerObject(1*AU, 0.01671022,
                     np.radians(0.00005),np.radians(-11.26064), np.radians(102.94719),
                     np.radians(100.46435), 5.97e24, m0, DTS*365.2)
Mercury = KeplerObject(0.38709893*AU, 0.20563069, np.radians(7.00487), np.radians(48.33167),
                       np.radians(77.45645), np.radians(252.25084), 0.33e24,m0, 88.0*DTS)
Venus = KeplerObject(0.72333199*AU, 0.00677323, np.radians(3.39471), np.radians(76.68069),
                     np.radians(131.53298), np.radians(181.97973),4.87e24, m0,224.7*DTS)
Mars = KeplerObject(1.52366231*AU, 0.09341233, np.radians(1.85061), np.radians(49.57854),
                     np.radians(336.04084), np.radians(355.45332),0.642e24, m0,687*DTS)
Jupiter = KeplerObject(5.20336301*AU, 0.04839266, np.radians(1.30530), np.radians(100.55615),
                       np.radians(14.75385), np.radians(34.40438),1898e24, m0,4331*DTS)
Saturn = KeplerObject(9.53707032*AU, 0.05415060, np.radians(2.48446), np.radians(113.71504),
                       np.radians(92.43194), np.radians(49.94432),568e24, m0,10747*DTS)
Uranus = KeplerObject(19.19126393*AU, 0.04716771, np.radians(0.76986), np.radians(74.22988),
                       np.radians(170.96424), np.radians(313.23218),86.8e24, m0,30589*DTS)
Neptune = KeplerObject(30.06896348*AU, 0.00858587, np.radians(1.76917), np.radians(131.72169),
                       np.radians(44.97135), np.radians(304.88003),102e24, m0,59800*DTS)
Pluto = KeplerObject(39.48168677*AU, 0.24880766, np.radians(17.14175), np.radians(110.30347),
                       np.radians(224.06676), np.radians(238.92881),0.0130e24, m0,90560*DTS)

solarsys = [Mercury, Venus, Earth, Mars,Jupiter,Saturn, Uranus, Neptune, Pluto]
#solarsys = [Earth]
print(Mercury.T / DT)



print("Hello!")
for i in range(N):  
    for obj in solarsys:
        obj.update(DT)
        data[i] += obj.pos.tolist()
        #print(obj.rotMatrix.as_matrix())    
   #break
   #data[i] = Earth.pos
   #print(Earth.pos)
#'''
data = np.transpose(np.array(data))
print(data.shape)

plotAll = True

#N = data.shape[1]
#print(data)
coords = [["X",0],["Y",1],["Z",2]]
toshow = [0,1]
numplanets = 2
colors = ["yellow","slategray","orange","deepskyblue", "firebrick",
              "darkgoldenrod","sandybrown","lightsteelblue","blue","dimgray"]
planets= ["Sun","Mercury","Venus","Earth","Mars",
              "Jupiter","Saturn","Uranus","Neptune","Pluto"]

numofplanets = len(solarsys)
#N = 100000
for i in range(0,N,10):
    break
    if plotAll:
        i = N
        K = N
    else:
        K = 40
    plt.figure(figsize= (6,6))
    for j in range(numofplanets):
        if (j > 4):
            continue
        plt.plot(data[3*j][max(0,i-K):i],data[3*j+1][max(0,i-K):i],color=colors[j+1], label = planets[j+1])
        
    plt.plot(0,0,'r+')
    plt.axis('square')
    #plt.xlim((-2e11,2e11))
    #plt.ylim((-2e11,2e11))
    #plt.plot(data[6][:i], data[7][:i], "y")
    plt.xlabel(coords[toshow[0]][0])
    plt.ylabel(coords[toshow[1]][0])
    plt.title("Positions in " + coords[toshow[0]][0] + coords[toshow[1]][0] + " plane for all time")
    plt.legend()
    plt.show()
    if plotAll:
        break
#'''

print(time.time() - starttime)