import matplotlib.pyplot as plt
from numpy import *
from random import choice
from copy import deepcopy
import sys
import time


if __name__=='__main__':


    directoryLoad = '/Users/Ren/Dropbox/SourceCode/test/MMPF/Result_inflow/'
    directoryLoad2 = '/Users/Ren/Dropbox/SourceCode/CORSIM filter factor/TrueState/'
    directoryLoad3 = '/Users/Ren/Dropbox/SourceCode/test/EnKF/Result_inflow/'
    directorySave = '/Users/Ren/Dropbox/SourceCode/test/MMPF/Result_PR_error/'

    TrueModel = load(directoryLoad+'TrueModel.npy')

    ErrorDensity = zeros((6,3))
    ErrorModel = zeros((6,3))

                
    currentRun = 1

    index=0
    for inflow in [1000,2000,3000,4000,5000,6000]:
        
        estDensity = load(directoryLoad + '4lag2estDen'+str(inflow)+'.npy')
        estModel = load(directoryLoad + '4lag2estModel'+str(inflow)+'.npy')
        estDensityPF = load(directoryLoad + '4lag0estDen'+str(inflow)+'PF.npy')
        estModelPF = load(directoryLoad + '4lag0estModel'+str(inflow)+'PF.npy')
        estDensityEnKF = load(directoryLoad3 + '4lag0estDen'+str(inflow)+'.npy')
        estModelEnKF = load(directoryLoad3 + '4lag0estModel'+str(inflow)+'.npy')

        TrueDensity = load(directoryLoad2 + 'TrueDensityL3B1F'+str(inflow)+'S65R50.npy')            
        TrueModel = load( directoryLoad + 'TrueModel.npy')

        ErrorDensity[index,0] = average(abs(estDensity - TrueDensity))
        ErrorModel[index,0] = average(abs(estModel - TrueModel))
        ErrorDensity[index,1] = average(abs(estDensityPF - TrueDensity))
        ErrorModel[index,1] = average(abs(estModelPF - TrueModel))
        ErrorDensity[index,2] = average(abs(estDensityEnKF - TrueDensity))
        ErrorModel[index,2] = average(abs(estModelEnKF - TrueModel))        
        index = index+1
        print 'flow is', inflow
        print 'MMPF',estModel[60:70,4]
        print 'EnKF',estModelEnKF[60:70,4]

    xData = [1000, 2000, 3000, 4000, 5000, 6000]

    print 'MMPF', ErrorDensity[:,0]
    print 'PF', ErrorDensity[:,1]
    print 'EnKF', ErrorDensity[:,2]

##    plt.rc('xtick',labelsize=20)
##    plt.rc('ytick',labelsize=20)
##    plt.hold(True)
##    StateNS=plt.plot(xData,ErrorDensity[:,0],marker='v',linestyle='-',color='b',label='$e_x$ MMPF')
##    StateWS=plt.plot(xData,ErrorDensity[:,1],marker='v',linestyle='-', color='r',label='$e_x$ PF')
##    StateEnKF=plt.plot(xData,ErrorDensity[:,2],marker='v',linestyle='-', color='g',label='$e_x$ IMM EnKF')
##    plt.legend(loc=1,prop={'size':20})
##    plt.xlim([0,7000])
##    plt.ylim([0,130])
##    plt.xticks([1000, 2000, 3000, 4000, 5000, 6000])
##    plt.xlabel('Inflow (veh/hour)',fontsize=20)
##    plt.ylabel('Error',fontsize=20)
##    plt.savefig(directorySave+'FigInflowError.pdf',bbox_inches='tight')
##    plt.show()
##    plt.hold(False)

##    plt.rc('xtick',labelsize=20)
##    plt.rc('ytick',labelsize=20)
##    plt.hold(True)
##    ParameterNS=plt.plot(xData,ErrorModel[:,0],marker='v',linestyle='--',color='b',label='$e_\gamma$ MMPF')
###    ParameterWS=plt.plot(xData,ErrorModel[:,1],marker='v',linestyle='--', color='r',label='$e_\gamma$ PF')
##    ParameterEnKF=plt.plot(xData,ErrorModel[:,2],marker='v',linestyle='--', color='g',label='$e_\gamma$ IMM EnKF')
##    plt.legend(loc=1,prop={'size':20})
##    plt.xlim([0,7000])
##    plt.ylim([0,0.05])
##    plt.xticks([1000, 2000, 3000, 4000, 5000, 6000])
##    plt.xlabel('Inflow (veh/hour)',fontsize=20)
##    plt.ylabel('Error',fontsize=20)
##    plt.savefig(directorySave+'FigInflowErrorModel.pdf',bbox_inches='tight')
##    plt.show()
##    plt.hold(False)

