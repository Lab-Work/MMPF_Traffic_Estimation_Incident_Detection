## This file calibrates the parameters for the fundamental diagram

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from numpy import *
from random import choice
from copy import deepcopy
import scipy.optimize as optimization
import sys
import os


if __name__=='__main__':

    ## the directory to load the processed keyData
    loadDirectory = os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir))+'/CORSIM_FD/ProcessedKeyData/'
    ## the directory to store the fundamental diagram
    saveDirectory = os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir))+'/CORSIM_FD/CalibratedFD/'       

    fileList = ['L3F5000S65','L3F3000S65','L3F7000S45','L3F7000S35','L3F7000S25','L3F7000S15','L3F7000S10','L3F7000S5','L3F7000S1']
    FDCollected = load(loadDirectory + 'L3F7000S65'+'ProcessedkeyData.npy')
    print FDCollected

    for fileName in fileList:
        FDCollectedFocus = load(loadDirectory + fileName+'ProcessedkeyData.npy')
        FDCollected = vstack((FDCollected,FDCollectedFocus))

    # data cleaning, remove observed noise points, the outlier point is decided by check the density--flow plot
    for k in range(len(FDCollected[:,0])):
            if 100<FDCollected[k,0]<150 and 700<FDCollected[k,2]<1200:
                print '=='
                print FDCollected[k,0]
                print FDCollected[k,2]
                FDCollected[k,0]=0.0
                FDCollected[k,2]=0.0
            if 20<FDCollected[k,0]<90 and 20<FDCollected[k,2]<800:
                print '=='
                print FDCollected[k,0]
                print FDCollected[k,2]
                FDCollected[k,0]=0.0
                FDCollected[k,2]=0.0


##################################### this section of codes use least square to fit the data ####################

    # Vmax is determined by CORSIM simulation
    Vmax = 65.0

    # Determine the capacity
    qmax = max(FDCollected[:,2])

    rhoc = qmax/Vmax

    def funC(x, a):
    # This function fits the data when the highway is congested
    # -b/(2a) = rhoc
        return a*x*x-2*rhoc*a*x+(4*a*qmax+4*a*a*rhoc*rhoc)/(4*a)

    # Separate the data to free flow data and congested data
    dataFreeFlow = zeros((2,sum(FDCollected[:,0]<rhoc)))
    dataCongested = zeros((2,sum(1-(FDCollected[:,0]<rhoc))))
    indexFreeFlow = 0
    indexCongested = 0

    for k in range(len(FDCollected[:,0])):
        if FDCollected[k,0]<rhoc:
            dataFreeFlow[0,indexFreeFlow] = FDCollected[k,0]
            dataFreeFlow[1,indexFreeFlow] = FDCollected[k,2]
            indexFreeFlow = indexFreeFlow+1 
        elif FDCollected[k,0]<260:
            dataCongested[0,indexCongested] = FDCollected[k,0]
            dataCongested[1,indexCongested] = FDCollected[k,2]
            indexCongested = indexCongested+1 

    # use least square fit to determine the parameters
    xc = [1]
    parameterCongested, errorCongested = optimization.curve_fit(funC, dataCongested[0,:], dataCongested[1,:], xc)
    a = parameterCongested[0]
    b = -2*rhoc*a
    c = (4*a*qmax+4*a*a*rhoc*rhoc)/(4*a)
        
    ## determine rhom
    for i in linspace(200,300,1000):
        j = i+(300-200)/1000.0
        if a*i*i+b*i+c >= 0 and a*j*j+b*j+c <= 0:
            rhom = i
            break

##################################### this section of codes generate the fundamental diagram ####################

    dataFlow = zeros((1000,2))
    dataSpeed = zeros((1000,2))

    n=0
    for i in linspace(0,250,1000):
        if i<rhoc:
            speed = Vmax
            flow = i*Vmax
        else:
            speed = a*i+b+c/i
            flow = a*i*i+b*i+c
        dataFlow[n] = [i,flow]
        dataSpeed[n] = [i,speed]
        n = n+1

    plt.rc('xtick',labelsize=30)
    plt.rc('ytick',labelsize=30)
    plt.hold(True)
    plt.scatter(FDCollected[:,0],FDCollected[:,1],s=4,color='r')
    plt.plot(dataFlow[:,0],dataSpeed[:,1],linewidth=2.5,color='k')
    plt.xlabel('Density (veh/mile/lane)',fontsize=30)
    plt.ylabel('Speed (veh/h/lane)',fontsize=30)
    plt.ylim([0,100])
    plt.xlim([0,300])
    plt.savefig(saveDirectory+'CORSIM_DS_'+str(int(rhoc))+'_'+str(int(rhom))+'.pdf',bbox_inches='tight')
    plt.show()
    plt.hold(False)

    plt.rc('xtick',labelsize=30)
    plt.rc('ytick',labelsize=30)
    plt.hold(True)
    plt.scatter(FDCollected[:,0],FDCollected[:,2],s=4,color='r')
    plt.plot(dataFlow[:,0],dataFlow[:,1],linewidth=2.5,color='k')
    plt.xlabel('Density (veh/mile/lane)',fontsize=30)
    plt.ylabel('Flow (veh/h/lane)',fontsize=30)
    plt.ylim([0,4000])
    plt.xlim([0,300])
    plt.savefig(saveDirectory+'CORSIM_DF_'+str(int(rhoc))+'_'+str(int(rhom))+'.pdf',bbox_inches='tight')
    plt.show()
    plt.hold(False)


