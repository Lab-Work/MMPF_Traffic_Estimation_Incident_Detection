import matplotlib.pyplot as plt
from numpy import *
from random import choice
from copy import deepcopy
import sys
import time


if __name__=='__main__':

    penetrationRange = [1,2,3,4]
    randomSample = 5
    smootherChoice = [True, False]
    savefigure = True

    directoryLoad = '/Users/Ren/Dropbox/SourceCode/test/MMPF/Result_all/'
    directorySave = '/Users/Ren/Dropbox/SourceCode/test/MMPF/Result_PR_error/'

    TrueDensity = load(directoryLoad+'TrueDensity.npy')
    TrueModel = load(directoryLoad+'TrueModel.npy')

    ErrorDensity = zeros((4,2))
    ErrorModel = zeros((4,2))

                
    currentRun = 1
    for PR in penetrationRange:
        for smoother in smootherChoice:
            for rand in range(randomSample):

#                print 'total 40, current:', currentRun
                currentRun = currentRun+1
                lag = 4
                if smoother == False:
                    lag=0


                estDensity = load( directoryLoad + str(PR)+'lag'+str(lag)+'R'+str(rand)+'estDen'+'.npy')
                estModel = load( directoryLoad + str(PR)+'lag'+str(lag)+'R'+str(rand)+'estModel'+'.npy')

                ErrorDensity[PR-1,smoother] = ErrorDensity[PR-1,smoother] + average(abs(estDensity - TrueDensity))
                ErrorModel[PR-1,smoother] = ErrorModel[PR-1,smoother] + average(abs(estModel - TrueModel))



    plt.rc('xtick',labelsize=20)
    plt.rc('ytick',labelsize=20)
    plt.hold(True)
    StateNS=plt.plot([1,2,3,4],ErrorDensity[:,0]/5,marker='v',linestyle='-',color='b',label='$e_x$ filter')
#    ParameterNS=plt.plot([1,2,3,4],ErrorModel[:,0]/5,marker='v',linestyle='--',color='b',label='$e_\gamma$ filter')
    StateWS=plt.plot([1,2,3,4],ErrorDensity[:,1]/5,marker='v',linestyle='--', color='r',label='$e_x$ smoother')
#    ParameterWS=plt.plot([1,2,3,4],ErrorModel[:,1]/5,marker='v',linestyle='--', color='r',label='$e_\gamma$ smoother')
    plt.legend(loc=1,prop={'size':20})
    plt.xlim([0.5,4.5])
    plt.ylim([0,50])
    plt.xticks([1.0,2.0,3.0,4.0])
    plt.xlabel('Penetration rate (%)',fontsize=20)
    plt.ylabel('Error',fontsize=20)
    plt.savefig(directorySave+'FigAllErrorDensity.pdf',bbox_inches='tight')
    plt.show()
    plt.hold(False)

    plt.rc('xtick',labelsize=20)
    plt.rc('ytick',labelsize=20)
    plt.hold(True)
#    StateNS=plt.plot([1,2,3,4],ErrorDensity[:,0]/5,marker='v',linestyle='-',color='b',label='$e_x$ filter')
    ParameterNS=plt.plot([1,2,3,4],ErrorModel[:,0]/5,marker='v',linestyle='-',color='b',label='$e_\gamma$ filter')
#    StateWS=plt.plot([1,2,3,4],ErrorDensity[:,1]/5,marker='v',linestyle='-', color='r',label='$e_x$ smoother')
    ParameterWS=plt.plot([1,2,3,4],ErrorModel[:,1]/5,marker='v',linestyle='--', color='r',label='$e_\gamma$ smoother')
    plt.legend(loc=1,prop={'size':20})
    plt.xlim([0.5,4.5])
    plt.ylim([0,0.05])
    plt.xticks([1.0,2.0,3.0,4.0])
    plt.xlabel('Penetration rate (%)',fontsize=20)
    plt.ylabel('Error',fontsize=20)
    plt.savefig(directorySave+'FigAllErrorModel.pdf',bbox_inches='tight')
    plt.show()
    plt.hold(False)







##                if PR == 2:
##                    print 'smoother', smoother
##                    print sum(abs(estModel - TrueModel))


                


##            ##########################################################################################################################
##            ## Load Measurements and true
##
##                directoryLoad = '/Users/Ren/Dropbox/SourceCode/CORSIM filter factor/'
##                directorySave = '/Users/Ren/Dropbox/SourceCode/test/MMPF/Result_occ/'
##
##
##                densityTrue = load(directoryLoad+'TrueState/'+'TrueDensityL3B1F6000S65R50.npy')
##            #    speedTrue = load(directoryLoad+'TrueSpeed.npy')
##                densityMea = densityTrue
##                speedMea = load(directoryLoad+'Measurements/'+'meaSpeed'+str(PR)+'R'+str(rand)+'L3B1F6000S65R50.npy')
##            #    parameterTrue = load(directoryLoad+'TrueParameter.npy')



