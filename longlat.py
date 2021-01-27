import numpy
import matplotlib.pyplot as plt
import sys
import os
import random
import math
import re
from scipy.optimize import linprog


def firstthings():
    '''Get data from files'''
    #latlong of census tracts
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
    
    #latlong of PODs
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

    #get population and vulnerable counts from census data
    tractpop = []
    tractvul = []
    file = open("/Users/kenyaandrews/Desktop/ResearchUIC/Fall2020/covid/censustract2020popvul.csv")
    for lines in file:
        lines = lines.replace('\t', ',').strip()
        lines = lines.replace(' ', ',')
        x, y = lines.split(',')
        tractpop.append(float(re.sub('[^0-9]','',x)))
        tractvul.append(float(y))
    file.close()
    print(tractpop, tractvul)

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
    alpha = 0.6
    #Number of vaccines for population tract
    N = [j * alpha for j in tractpop]

    #calculate Vaccination Disparity Vulnerable and Non-Vulnerable Rate => star
    #Calculate Hospitalization Rate Vul vs Non-Vul => doublestar
    star = [] #Vaccination Disparity Vulnerable and Non-Vulnerable Rate
    nonvulrate = []
    doublestar = []
    rho = numpy.zeros((len(eta), len(tractpop)))
    #linev = [0.6-0.6*j for j in eta]
    gamma = 0
    n_i= c_i = 0
    c_iarray = numpy.zeros((len(tractpop),1))
    epsilon = 0.1
    l1_bound = 0
    totalpop = sum(tractpop)
    
    totalVaccines = alpha*totalpop
    for i in range(len(eta)):
        vulnerablestar = nonvulnerablestar = 0.0
        numbervul = 0.0
        numbernonvul = 0.0
        for j in range(len(tractpop)):
            x=(percentvul[j]*eta[i])/(percentvul[j]*eta[i]+1-percentvul[j]) #rho_i
            rho[i][j] = x
            vulnerablestar += x*N[j] #rho_i * N_i
            numbervul += percentvul[j]*tractpop[j] #beta_i * C_i
            
            y=(1-percentvul[j])/(percentvul[j]*eta[i]+1-percentvul[j])
            nonvulnerablestar += y*N[j]
            
            numbernonvul += (1-percentvul[j])*tractpop[j]
            #print((x+y))
            gamma += numbervul/totalpop
            n_i= (N[j]/alpha)/totalVaccines
            c_i= tractpop[j]/totalpop
            c_iarray[j] = c_i
            #minimize CNvrd
            l1_bound += abs(n_i- c_i)
            CNvrd = ((-1/gamma)*(-vulnerablestar/numbervul)*n_i) + (((1-(1-gamma))* nonvulnerablestar/numbernonvul)*n_i)
        
        ehospitalvul = ((5*(numbervul-vulnerablestar)) + (0.5*vulnerablestar))
        ehospitalnonvul = ((numbernonvul-nonvulnerablestar) + (0.1*nonvulnerablestar))
        
        nonvulrate.append(nonvulnerablestar/numbernonvul)
        
        
        star.append(-vulnerablestar/numbervul + nonvulnerablestar/numbernonvul)
        doublestar.append((ehospitalvul/numbervul)/(ehospitalnonvul/numbernonvul))
    
    '''TESTING HERE'''
    #print(star[0])
    #print("Lower Bound", l1_bound)
    #print(eta)
    #print(nonvulrate[0])
    #print(len([N[j]/tractpop[j] for j in range(len(tractpop))]))
    #print(percentvul)
    #print(numbernonvul/(numbernonvul+numbervul))
    #print(0.6/(numbernonvul/(numbernonvul+numbervul)))
    '''END TESTS'''
    
    #plot Varying Eta vs Vaccination Disparity Vulnerable and Non-Vulnerable Rate
    fig1, starplot = plt.subplots()
    starplot.plot(eta, star)
    plt.xlabel('Eta')
    plt.ylabel('Disparity Value (Vulnerable vs Non-Vulnerable)')
    plt.title('Varying Eta vs Vaccination Disparity Vulnerable and Non-Vulnerable Rate')
    #plt.show()
    
    fig2, doublestarplot = plt.subplots()
    doublestarplot.plot(eta, star)
    plt.xlabel('Eta')
    plt.ylabel('Expected Hospitalization Rate Disparity')
    plt.title('Varying Eta vs Expected Hospitalization Rate Disparity')
    #plt.show()
  
    k = len(tractpop)
    #get A
    A = numpy.zeros( (3*k+1, 2*k) )
    
    counter = 0
    row = len(A)
    for i in range(k):
        A[i][counter] = 1
        A[i][counter+k] = -1
        
        A[i+k][counter] = -1
        A[i+k][counter+k] = -1
        
        A[i+(2*k)][counter+k] = -1
     
        counter+=1

    for i in range(int(round(A.shape[1]/2)),):
        A[A.shape[0]-1][i+int(round(A.shape[1]/2))] = 1
    
    
    #get b
    b = numpy.zeros( (3*k+1, 1) )
    b[:k] = c_iarray
    b[k:2*k] = -c_iarray
    b[-1] = epsilon
    
    #get Aeq
    Aeq = numpy.zeros( (1, 2*k) )
    
    for i in range(int(round(Aeq.shape[1]/2)),):
        Aeq[0][i+int(round(Aeq.shape[0]/2))] = 1
    print("Aeq", Aeq)
    c_ = numpy.zeros((2*k, 1))
    beq = numpy.zeros((1, 1))
    #set beq
    beq[0] = 1

    
    #try the optimzation with scipy.linprog
    print("Optimzation: ")

    #print("Size of c: ", c_.shape)
    print("Size of A: ", A.shape)
    print("Size of b: ", b.shape)
    print("Size of beq: ", beq.shape)
    print("Size of Aeq: ", Aeq.shape)
    
    #get c
    c = (-(1/gamma)*rho[5]) + ((1/(1-gamma))*(1-rho[5]))
    c_ = numpy.zeros((2*k, 1))
    c_[:k] = c.reshape(k,1)
    
    res = linprog(c = c_, A_ub=A, b_ub=b, A_eq = Aeq, b_eq = beq, options={"disp": True})
    #res = linprog(c = numpy.array(c_), A_ub=numpy.array(A), b_ub=numpy.array(b),A_eq = Aeq, b_eq = beq, options={"disp": True})
    #print(res)
    print(res.x)
    print("sum resx", sum(res.x[:k]))
    x = numpy.zeros((2*k, 1))
    x[:k] = c_iarray
    '''
    print("A*x", numpy.matmul(A,x)-b )
    print("aeq*x", numpy.matmul(Aeq,x)-beq)
    print("c_is", numpy.sum(c_iarray))
    '''
    
    #Number of vaccines for population tract
     #N = [j * 0.6 for j in tractpop]

     #calculate Vaccination Disparity Vulnerable and Non-Vulnerable Rate => star
     #Calculate Hospitalization Rate Vul vs Non-Vul => doublestar
    star = [] #Vaccination Disparity Vulnerable and Non-Vulnerable Rate
    nonvulrate = []
    doublestar = []
    epsilon = 0.1
    for i in range(len(eta)):
         vulnerablestar = nonvulnerablestar = 0.0
         numbervul = 0.0
         numbernonvul = 0.0
         #get c
         c = (-(1/gamma)*rho[i]) + ((1/(1-gamma))*(1-rho[i]))
         c_ = numpy.zeros((2*k, 1))
         c_[:k] = c.reshape(k,1)
         
         res = linprog(c = c_, A_ub=A, b_ub=b, A_eq = Aeq, b_eq = beq, options={"disp": True})
         N = res.x[:k]*totalVaccines
         print("any less :", any(N<0))
         for j in range(len(tractpop)):
             x=(percentvul[j]*eta[i])/(percentvul[j]*eta[i]+1-percentvul[j]) #rho_i
             
             vulnerablestar += x*N[j] #rho_i * N_i
             numbervul += percentvul[j]*tractpop[j] #beta_i * C_i
             
             y=(1-percentvul[j])/(percentvul[j]*eta[i]+1-percentvul[j])
             nonvulnerablestar += y*N[j]
             
             numbernonvul += (1-percentvul[j])*tractpop[j]
         
         ehospitalvul = ((5*(numbervul-vulnerablestar)) + (0.5*vulnerablestar))
         ehospitalnonvul = ((numbernonvul-nonvulnerablestar) + (0.1*nonvulnerablestar))
         
         nonvulrate.append(nonvulnerablestar/numbernonvul)
         
         
         star.append(-vulnerablestar/numbervul + nonvulnerablestar/numbernonvul)
         doublestar.append((ehospitalvul/numbervul)/(ehospitalnonvul/numbernonvul))


    #plot Varying Eta vs Vaccination Disparity Vulnerable and Non-Vulnerable Rate
    starplot.plot(eta, star)
    plt.show()

    doublestarplot.plot(eta, doublestar)
    plt.show()
firstthings()









