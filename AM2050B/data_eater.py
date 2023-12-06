
import matplotlib.pyplot as plt
import numpy as np
import time
import os

datafiles =[ os.path.join('data/',x) for x in os.listdir('data/')]
newestfile = max(datafiles, key = os.path.getctime)


with open(newestfile) as f:
    rawdata = f.readlines()

data = [[float(j) for j in rawdata[i].replace(" ", "").replace("\n","").split(",")] for i in range(len(rawdata))]
numplanets = len(data[0]) // 3
data = np.transpose(np.array(data))
print(data)

T = 20  # seconds
N = 20000
DT = T/N
funny = False
for i in range(9000,N,10): 
    #'''
    if (not funny and i > N/2):
        funny = True
        print("WE'RE GOING BACK! BACK TO THE FUTURE!")
    #'''
    #i = N
    K = 100
    plt.figure(figsize= (6,6))
    colors = ["mediumseagreen","cornflowerblue", "yellow","salmon","hotpink", "deepskyblue","orange"]
    for j in range(numplanets):
        if (j >= len(colors) or j >= 2):
            break
        plt.plot(data[3*j][max(0,i-K):i], data[3*j+1][max(0,i-K):i], color=colors[j])
        
    #plt.plot(data[3][0],data[4][0],'r+')
    #plt.plot(data[6][:i], data[7][:i], "y")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Positions in XY plane for all time")
    plt.show()
    #break

    #time.sleep(DT)
