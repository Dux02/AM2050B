
import matplotlib.pyplot as plt
import numpy as np
import time
import os

datafiles =[ os.path.join('data/',x) for x in os.listdir('data/')]
newestfile = max(datafiles, key = os.path.getctime)

if (newestfile.startswith("data/kss")):
    print("Yipee")
    colors = ["yellow","slategray","orange","deepskyblue", "firebrick",
              "darkgoldenrod","sandybrown","lightsteelblue","blue","dimgray"]
    planets= ["Sun","Mercury","Venus","Earth","Mars",
              "Jupiter","Saturn","Uranus","Neptune","Pluto"]
else:
    colors = ["mediumseagreen","cornflowerblue", "yellow","salmon","hotpink", "deepskyblue","orange"]
    planets = ["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn"]
    
plotAll = True

with open(newestfile) as f:
    rawdata = f.readlines()

data = [[float(j) for j in rawdata[i].replace(" ", "").replace("\n","").split(",")] for i in range(len(rawdata))]
numplanets = len(data[0]) // 3
data = np.transpose(np.array(data))
N = data.shape[1]
#print(data)
coords = [["X",0],["Y",1],["Z",2]]
toshow = [0,1]

#N = 100000
for i in range(0,N,100): 
    if plotAll:
        i = N
        K = N
    else:
        K = 400
    plt.figure(figsize= (6,6))
    for j in range(numplanets):
        if (j >= len(colors) or j == 0 or j > 6):
            continue
        plt.plot(data[3*j + coords[toshow[0]][1]][max(0,i-K):i], 
                 data[3*j+ coords[toshow[1]][1]][max(0,i-K):i], 
                 color=colors[j],label=planets[j])
        
    plt.plot(data[0][0],data[1][0],'r+')
    plt.axis('square')
    #plt.xlim((-1e12,1e12))
    #plt.ylim((-1e12,1e12))
    #plt.plot(data[6][:i], data[7][:i], "y")
    plt.xlabel(coords[toshow[0]][0])
    plt.ylabel(coords[toshow[1]][0])
    plt.title("Positions in " + coords[toshow[0]][0] + coords[toshow[1]][0] + " plane for all time")
    plt.legend()
    plt.show()
    if plotAll:
        break

    #time.sleep(DT)
