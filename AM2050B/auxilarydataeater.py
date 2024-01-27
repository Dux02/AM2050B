# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 00:27:23 2024

@author: Blas
"""

rawdata = open("data/dtit 01-27 00-22.txt").readlines()
import matplotlib.pyplot as plt
import numpy as np

N = len(rawdata)
data = [[float(j) for j in rawdata[i].replace(" ", "").replace("\n","").split(",")] for i in range(int(0*N),N,1)]
data = np.transpose(np.array(data))
#T, Earth xyz, Asteroid xyz × 3

def plotplanetstuff():

    plt.figure(figsize =(10,6))
    def norm(j):
        return np.sqrt(data[3*j+4]*data[3*j+4] + data[3*j+5]*data[3*j+5] + 
                       data[3*j+6]*data[3*j+6])
    plt.plot(data[0],norm(0),'k.')
    plt.plot(data[0],norm(1),'r.')
    plt.plot(data[0],norm(2),'b.')
    print(data[0][-1],norm(0)[-1],norm(1)[-1], norm(2)[-1])
    print("Position DT=640",data[4][-1],data[5][-1],data[6][-1])
    print("Position DT=3600",data[7][-1],data[8][-1],data[9][-1])
    print("Position DT=2880",data[10][-1],data[11][-1],data[12][-1])
    plt.yscale("log")
    plt.xlabel("T (s)")
    plt.ylabel("Distance from sun (m)")
    plt.show()
    colors = ["deepskyblue","k","r","b"]
    planets = ["Earth","DT = 640","DT = 3600", "DT = 2880"]
    
    plotAll = True
    DT = 365.2*24*60*60
    
    for i in range(0,N,100):
        if plotAll:
            i = N
            K = N
        else:
            K = 400
        
        plt.figure(figsize = (6,6))
        for j in range(4):
            if (plotAll):
                pass
            plt.plot(data[3*j + 1][max(0,i-K):i], 
                     data[3*j + 2][max(0,i-K):i], 
                     color=colors[j],label=planets[j])
            plt.plot(data[3*j + 1][i-1], 
                     data[3*j + 2][i-1],'+',
                     color=colors[j],ms=10)
        plt.plot(1.8924675e10, -1.4591901e11,'r.',ms=10)
        plt.axis('square')
        titlestr = "Positions in XY plane"
        if plotAll:
            titlestr += " for all time T= " + str(round(data[0][N//2]/DT)) + " years"
        else:
            titlestr += " at time T= -" + str(round(data[0][i]/DT,2)) + " years"
        plt.title(titlestr)
        plt.legend()
        plt.show()
        if plotAll:
            break
def plottimestuff():
    plt.figure(figsize =(10,6))
    plt.plot(data[0],data[1],'k.')
    plt.yscale("log")
    plt.xlabel("DT (s)")
    plt.ylabel("Absolute difference (m)")
    plt.show()
    print(data[0][np.argmax(data[1])])

plotplanetstuff()