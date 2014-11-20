import matplotlib.pyplot as plt
from numpy import *
from random import choice
from copy import deepcopy
import sys
sys.path.append('/Users/Ren/Dropbox/SourceCode/FilterFunctions/')
sys.path.append('/Users/Ren/Dropbox/SourceCode/TrafficModel/')
from traffic_models_lq import *
from kalman_filter_functions import *
import time

####################################  statictical boundary conditions   #######################################

def left_boundary_statistic(laneNumber, Vmax):
    lbDensity = (6900.0+ random.normal(0,100))/Vmax
    lbModel = laneNumber
    return lbDensity, lbModel

def right_boundary_statistic(laneNumber, Vmax):
    rbDensity = (6900.0+ random.normal(0,100))/Vmax
    rbModel = laneNumber
    return rbDensity, rbModel

def init_state(DensityMeaInit,cellNumber,laneNumber, sample):
    mean = (DensityMeaInit[1]+DensityMeaInit[-2])/2
    std = 0.05*mean
    state = random.normal(mean, std, (sample,cellNumber))
    return state


###############################################  plot functions   ##################################################


def plot_density(data, bool, savefile, directorySave):
    plt.rc('xtick',labelsize=30)
    plt.rc('ytick',labelsize=30)
    plt.imshow(data,aspect='auto',origin='lower',interpolation='nearest')
    plt.ylabel('Time Step',fontsize=30)
    plt.clim(0.0, 560)
    plt.xlabel('Cell Number',fontsize=30)
    plt.colorbar()
    if bool == True:
        plt.savefig(directorySave + savefile+'.pdf', bbox_inches='tight')
    plt.show()
    plt.clf()

def plot_parameter(data, bool, savefile, directorySave):
    cmap=plt.cm.jet_r
    plt.rc('xtick',labelsize=30)
    plt.rc('ytick',labelsize=30)
    plt.imshow(data,aspect='auto',origin='lower',interpolation='nearest',cmap=cmap)
    plt.clim(0,3)
    plt.ylabel('Time Step',fontsize=30)
    plt.xlabel('Cell Number',fontsize=30)
    plt.colorbar()
    if bool == True:
        plt.savefig(directorySave + savefile+'.pdf', bbox_inches='tight')
    plt.show()
    plt.clf()
    
##########################################################################################################################

if __name__=='__main__':

##  Model noise and measurement noise
    modelNoiseMean = 0.0
    modelNoiseStd = 5.0

    densityMeaMean = 0.0
    densityMeaStd = 13.5 #13.5

    speedMeaMean = -4.0
    speedMeaStd = 4.8

##  Parameters for traffic models

    rhoc = 34.0
    rhom = 239.0    
    Vmax = 65.0
#    rhom_tuide = 1303.0

    rhocOneIncident = 90.2 #115.0
    rhocTwoIncident = 62.6
    VmaxIncident = 18.0
#    rhom_tuideIncident = 30000.0


##  Parameters for simulation
    dt = 20.0 / 3600
    length = 4.0
    laneNumber = 3
    maxLaneBlocked = 3
    sample = 50
    timeStep = 180
    pTran = 0.99

    PR = 4
    savefigure = True

##  Discretization
    dx = Vmax*dt
    cellNumber = floor(length/dx)
    dx = length/cellNumber
    Lambda = dt/dx

## Load Measurements and true

    directoryLoad = '/Users/Ren/Dropbox/SourceCode/CORSIM filter factor/'
    directorySave = '/Users/Ren/Dropbox/SourceCode/test/EnKF/result/'


    densityTrue = load(directoryLoad+'TrueState/'+'TrueDensityL3B1F6000S65R50.npy')
#    speedTrue = load(directoryLoad+'TrueSpeed.npy')
    densityMea = densityTrue
    speedMea = load(directoryLoad+'Measurements/'+'meaSpeed'+str(PR)+'R0L3B1F6000S65R50.npy')
    parameterTrue = load(directoryLoad+'TrueState/'+'TrueModel.npy')


## compuate mode, model number
    modelSet, modelNumber = build_system_models_nb(cellNumber, laneNumber, maxLaneBlocked)

## Creat array to save results
    estimatedState = zeros((timeStep, cellNumber))
    estimatedModel = zeros((timeStep, cellNumber))
    
    estimatedStateAllModel = zeros((modelNumber, sample, cellNumber))
    modelProbability = zeros(modelNumber)
    modelLikelihood = zeros(modelNumber)

## Initialization

    state = init_state(densityMea, cellNumber, laneNumber, sample)
    speed = Vmax*ones((sample,cellNumber))

    estimatedState[0] = average(state, axis=0)
    estimatedModel[0] = laneNumber*ones(cellNumber)

##########################################################################################################################

## Algorithm Start
    start_time = time.time()

    for counter in range(1,timeStep):
        print 'this is interation', counter

# initialize model probability
        previousModel = estimatedModel[counter-1]
        modelProbability = init_model_probability_nb(previousModel, pTran, modelNumber, laneNumber, maxLaneBlocked)
        stateCopy = state.copy()

        densityMeaStep = zeros(cellNumber)
        densityMeaStep[1] = densityMea[counter, 1]
        densityMeaStep[-2] = densityMea[counter, -2]
        speedMeaStep = speedMea[counter]
        
        boundaryLeft = left_boundary_statistic(laneNumber, Vmax)
        boundaryRight = right_boundary_statistic(laneNumber, Vmax)
        


##        for modelIndex in range(0,1):
##            model = laneNumber*ones(cellNumber)
        for modelIndex in range(0,int(modelNumber)):
            model = modelSet[modelIndex]
            state = stateCopy.copy() 
# forward prediction
            for j in range(int(sample)):
                state[j] = ctm_lq(state[j],model, Lambda, boundaryLeft, boundaryRight, rhoc, rhom, Vmax, modelNoiseMean, modelNoiseStd, rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber)
                speed[j] = vel_lq(state[j], rhoc, rhom, Vmax, model, rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber)
                
# kalman gain
            kalmanGain = compute_KalmanGain(state, speed, cellNumber, sample, densityMeaStep, speedMeaStep, densityMeaStd, speedMeaStd)   

# state estimation
            state = kalman_state_update(state, speed, cellNumber, sample, densityMeaStep, speedMeaStep, kalmanGain, densityMeaMean, densityMeaStd, speedMeaMean, speedMeaStd)

# model likelihood
            modelLikelihood[modelIndex] = compute_model_likelihood(state, speed, cellNumber, densityMeaStep, speedMeaStep, densityMeaMean, densityMeaStd, speedMeaMean, speedMeaStd)
            

            estimatedStateAllModel[modelIndex] = state.copy()
            
### update model probability
        modelLikelihood = modelLikelihood/sum(modelLikelihood)
        print modelLikelihood
        modelProbability = modelProbability*modelLikelihood
        print modelProbability
        optimalIndex = modelProbability.argmax()
        state = estimatedStateAllModel[optimalIndex]
        estimatedModel[counter]=modelSet[optimalIndex]
        estimatedState[counter]=average(state, axis=0)
        print estimatedModel[counter]

##########################################################################################################################

    stateError = estimatedState - densityTrue[0:180]
    modelError = estimatedModel - parameterTrue[0:180]

    print "run time is", time.time() - start_time
    print 'density error is', average(abs(stateError))
    print 'model error is', sum(abs(modelError))

##    save( directorySave + str(PR)+'estDen',estimatedState)
##    save( directorySave + str(PR)+'estModel',estimatedModel)
##
##    plot_density(densityTrue, savefigure, 'trueDensityPlot', directorySave)
##    plot_density(estimatedState, savefigure, str(PR)+'estDen', directorySave)
##    plot_parameter(estimatedModel, savefigure, str(PR)+'EstModel', directorySave)
#    plot_parameter(parameterTrue[0:180], savefigure, 'trueModelPlot', directorySave)







