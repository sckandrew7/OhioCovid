import numpy
import matplotlib.pyplot as plt
import sys
import os
import random
import math

#PRINT THE PLOT WITHOUT THE CENTROIDS
def firstthings():
    file = open("/Users/kenyaandrews/Desktop/ResearchUIC/Fall2020/covid/censuslatlong.csv") #open A.txt
    censustractx = []
    censustracty = []
    for lines in file:
        lines = lines.replace('\t', ',').strip()
        lines = lines.replace(' ', ',')
        x, y = lines.split(',')
        censustractx.append(float(x))
        censustracty.append(float(y))
    file.close()
    theogplot = plt.scatter(censustractx, censustracty, color="orange")
    file = open("/Users/kenyaandrews/Desktop/ResearchUIC/Fall2020/covid/podlatlong.csv") #open A.txt
    podx = []
    pody = []
    for lines in file:
        lines = lines.replace('\t', ',').strip()
        lines = lines.replace(' ', ',')
        x, y = lines.split(',')
        podx.append(float(x))
        pody.append(float(y))
    file.close()

    theogplot = plt.scatter(podx, pody)
    plt.show()

#get pop and vul data
    tractpop = []
    tractvul = []
    file = open("/Users/kenyaandrews/Desktop/ResearchUIC/Fall2020/covid/censustract2010popvul.csv")
    for lines in file:
        lines = lines.replace('\t', ',').strip()
        lines = lines.replace(' ', ',')
        x, y = lines.split(',')
        tractpop.append(float(x))
        tractvul.append(float(y))
    file.close()
    
    #ASSIGN THE CLOSEST POINTS
    closestpodcensusx = []
    closestpodcensusy = []
    podpopsum = numpy.zeros(len(podx))
    podvulsum = numpy.zeros(len(podx))
    dismatrix = numpy.empty(shape=(len(censustractx),len(podx)))

    podpopsum[0] = tractpop[0]


    for i in range(len(censustractx)):
        for j in range(len(podx)):
            dismatrix[i][j] = math.sqrt((podx[j]-censustractx[i])**2 + (pody[j]-censustracty[i])**2)
        cen, podsmall = numpy.unravel_index(dismatrix[i].argmin(), dismatrix.shape)
        closestpodcensusx.append(podx[podsmall])
        closestpodcensusy.append(pody[podsmall])
        print("smallest for (", censustractx[i] , censustracty[i], ")  (", podx[podsmall], pody[podsmall], ")")
    k=0
    for i in range(1, len(podx)-1):
        if closestpodcensusx[i-1] == closestpodcensusx[i] and closestpodcensusy[i-1] == closestpodcensusy[i]:
            podpopsum[k] += tractpop[i]
            podvulsum[k] += tractvul[i]
        else:
            k+=1
            podpopsum[k] = tractpop[i]
            podvulsum[k] = tractvul[i]

    print(podpopsum)
    print(podvulsum)
        #we want to get number vulnerable and the population size for each pod
    #take a pod, and the censustract that goes with it. find the population/vulnerable and add it


firstthings()









