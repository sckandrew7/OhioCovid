import numpy
import matplotlib.pyplot as plt
import sys
import os
import random
import math


def firstthings():
    '''Get data from files'''
    file = open("/Users/kenyaandrews/Desktop/ResearchUIC/Fall2020/covid/censuslatlong.csv") #open census lat/long only
    censustractlat = []
    censustractlong = []
    for lines in file:
        lines = lines.replace('\t', ',').strip()
        lines = lines.replace(' ', ',')
        x, y = lines.split(',')
        censustractlat.append(float(x))
        censustractlong.append(float(y))
    file.close()
    theogplot = plt.scatter(censustractlong, censustractlat, color="orange")
    file = open("/Users/kenyaandrews/Desktop/ResearchUIC/Fall2020/covid/podlatlong.csv") #open pods lat/long only
    podlat = []
    podlong = []
    for lines in file:
        lines = lines.replace('\t', ',').strip()
        lines = lines.replace(' ', ',')
        x, y = lines.split(',')
        podlat.append(float(x))
        podlong.append(float(y))
    file.close()

    theogplot = plt.scatter(podlong, podlat)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Census Tract Centers (Orange) and PODs (Blue)')
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
    
    '''ASSIGN THE CLOSEST PODs'''
    closestpodcensuslat = []
    closestpodcensuslong = []
    podpopsum = numpy.zeros(len(podlat))
    podvulsum = numpy.zeros(len(podlat))
    dismatrix = numpy.empty(shape=(len(censustractlat),len(podlat)))


    for i in range(len(censustractlat)):
        for j in range(len(podlat)):
            dismatrix[i][j] = math.sqrt((podlat[j]-censustractlat[i])**2 + (podlong[j]-censustractlong[i])**2)
        cen, podsmall = numpy.unravel_index(dismatrix[i].argmin(), dismatrix.shape)
        closestpodcensuslat.append(podlat[podsmall])
        closestpodcensuslong.append(podlong[podsmall])
        print("smallest for (", censustractlat[i] , censustractlong[i], ")  (", podlat[podsmall], podlong[podsmall], ")")


    '''print the populations size and vulnerable size to their tracts, only storing values with a POD near them'''
    podpopsum[0] = tractpop[0]
    trackoffpointslat = []
    trackoffpointslong = []
    trackoffpointslat.append(closestpodcensuslat[0])
    trackoffpointslong.append(closestpodcensuslong[0])
    k=0
    for i in range(1, len(podlat)-1):
        if closestpodcensuslat[i-1] == closestpodcensuslat[i] and closestpodcensuslong[i-1] == closestpodcensuslong[i]:
            podpopsum[trackoffpointslat.index(closestpodcensuslat[i])] += tractpop[i]
            podvulsum[trackoffpointslat.index(closestpodcensuslat[i])] += tractvul[i]
        else:
            k+=1
            trackoffpointslat.append(closestpodcensuslat[i])
            trackoffpointslong.append(closestpodcensuslong[i])
            podpopsum[k] = tractpop[i]
            podvulsum[k] = tractvul[i]

    #only print population and vulnerable count for PODs closest to tract centers
    for i in range(len(trackoffpointslat)):
        print("Total Population for", trackoffpointslong[i],trackoffpointslat[i], "-" , podpopsum[i])
        print("Vulnerable Population for", trackoffpointslong[i],trackoffpointslat[i], "-" , podvulsum[i])


    '''graphing Varying Eta vs Vaccination Disparity Vulnerable and Non-Vulnerable Rate'''

    #bias
    eta = numpy.arange(0.0, 1.1, 0.1)
    
    #calculate percent vul
    percentvul = [] #Beta
    for i in range(len(tractpop)):
        percentvul.append(tractvul[i]/tractpop[i]) #CALCULATE BETA(percentage vulnerable)

    #Number of vaccines for population tract
    N = [j * 0.6 for j in tractpop]
    
    #calculate Vaccination Disparity Vulnerable and Non-Vulnerable Rate = star
    star = [] #Vaccination Disparity Vulnerable and Non-Vulnerable Rate
    vulnerablestar = nonvulnerablestar = 0
    for i in range(len(eta)):
        for j in range(len(tractpop)):
            vulnerablestar += ((percentvul[j]*eta[i]*N[j])/(percentvul[j]*eta[i]+1-percentvul[j]))/(percentvul[j]*tractpop[j])
            nonvulnerablestar += (((1-percentvul[j])*N[j])/(percentvul[j]*eta[i]+1-percentvul[j]))/((1-percentvul[j])*tractpop[j])
        star.append(-vulnerablestar + nonvulnerablestar)
    
    #plot Varying Eta vs Vaccination Disparity Vulnerable and Non-Vulnerable Rate
    starplot = plt.plot(eta, star)
    plt.xlabel('Eta')
    plt.ylabel('Disparity Value (Vulnerable vs Non-Vulnerable)')
    plt.title('Varying Eta vs Vaccination Disparity Vulnerable and Non-Vulnerable Rate')
    plt.show()
  


firstthings()









