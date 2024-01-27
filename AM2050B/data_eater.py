
import matplotlib.pyplot as plt
import numpy as np
import time
import os

datafiles =[ os.path.join('data/',x) for x in os.listdir('data/')]
newestfile = max(datafiles, key = os.path.getctime)

if newestfile.startswith("data/kss"):
    print("Yipee")
    colors = ["yellow","slategray","orange","deepskyblue", "firebrick",
              "darkgoldenrod","sandybrown","lightsteelblue","blue","dimgray"]
    planets= ["Sun","Mercury","Venus","Earth","Mars",
              "Jupiter","Saturn","Uranus","Neptune","Pluto"]
elif newestfile == "data/collision backy.txt" or newestfile.startswith('data/col'):
    colors = ["yellow","red","slategray","orange","deepskyblue", "firebrick",
              "darkgoldenrod","sandybrown","lightsteelblue","blue","dimgray"]
    planets= ["Sun","Asteroid","Mercury","Venus","Earth","Mars",
              "Jupiter","Saturn","Uranus","Neptune","Pluto"]
    
else:
    colors = ["mediumseagreen","cornflowerblue", "yellow","salmon","hotpink", "deepskyblue","orange"]
    planets = ["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn"]
    
plotAll = False
print(newestfile)

with open(newestfile) as f:
    rawdata = f.readlines()

N = len(rawdata)
data = [[float(j) for j in rawdata[i].replace(" ", "").replace("\n","").split(",")] for i in range(int(0*N),N,1)]
numplanets = len(data[0]) // 3
data = np.transpose(np.array(data))
N = int(data.shape[1])
#print(data)
coords = [["X",0],["Y",1],["Z",2]]
toshow = [0,1]
DT = 640 / (60*60*24)

indexsplindex = [0,10,20,30,40,50,60,70,80,90,100]
#N = 100000
rewound = False
for i in range(int(0/DT),N,1):
    if (not (i in indexsplindex or N-i in indexsplindex)):
        continue
    print("Booyah!")
    if (i >= int(N/2) and not rewound):
        print("Rewinding time")
        rewound = True
    #break
    if plotAll:
        i = N
        K = N
    else:
        K = 4000
    plt.figure(figsize= (6,6))
    for j in range(numplanets):
        if (j >= len(colors) or j == 0 or not (j in [1,2,3,4])):
            continue
        if (plotAll):
            pass
        elif (j == 1):
            K = 4000
            pass
        else:
            K = 4000
            pass
        plt.plot(data[3*j + coords[toshow[0]][1]][max(0,i-K):i], 
                 data[3*j+ coords[toshow[1]][1]][max(0,i-K):i], 
                 color=colors[j],label=planets[j])
        plt.plot(data[3*j + coords[toshow[0]][1]][i-1], 
                 data[3*j+ coords[toshow[1]][1]][i-1],'+',
                 color=colors[j],ms=10)
        
    plt.plot(data[0][0],data[1][0],'r+')
    plt.plot(1.8924675e10, -1.4591901e11,'r.',ms=10)
    plt.axis('square')
    plt.xlim((0.184e11,0.194e11))
    plt.ylim((-1.45e11,-1.46e11))
    #plt.plot(data[6][:i], data[7][:i], "y")
    plt.xlabel(coords[toshow[0]][0])
    plt.ylabel(coords[toshow[1]][0])
    titlestr = "Positions in " + coords[toshow[0]][0] + coords[toshow[1]][0] + " plane"
    if plotAll:
        titlestr += " for all time T= -" + str(round(DT*N)) + " days"
    else:
        if (rewound):
            titlestr += " at time T= -" + str(round(DT*(N-i),2)) + " days (rewound)"
        else:
            titlestr += " at time T= -" + str(round(DT*i,2)) + " days"
    plt.title(titlestr)
    plt.legend()
    plt.show()
    if plotAll:
        break

    #time.sleep(DT)
'''
rawdata = open('rotmat.txt').readlines()
interest = [[float(j) for j in rawdata[i].replace(" ", "").replace("\n","").split(",")] for i in range(len(rawdata))]
interest = np.transpose(np.array(interest))
plt.plot(interest[0],interest[1],'k-')
print(interest[1][20000])
print(interest[1][19999])
print(interest[1][0])
plt.axis('square')
#E = np.linspace(np.min(interest[1]), np.max(interest[1]),100)
#plt.plot(E-0.04839266*np.sin(E),E,'r--')
plt.show()
#'''