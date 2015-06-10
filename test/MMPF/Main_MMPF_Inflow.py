import matplotlib.pyplot as plt
from numpy import *
from random import choice
from copy import deepcopy
import sys
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir))+'/FilterFunctions/')
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir))+'/TrafficModel/')
from traffic_models_lq import *
from particle_filter_functions import *
import time

####################################  statictical boundary conditions   #######################################

def left_boundary_statistic(laneNumber, Vmax, inflow):
    lbDensity = (inflow-100 + random.normal(0,150))/Vmax
    lbModel = laneNumber
    return lbDensity, lbModel

def right_boundary_statistic(laneNumber, Vmax):
    rbDensity = (5900.0+ random.normal(0,150))/Vmax
    rbModel = laneNumber
    return rbDensity, rbModel

def init_state(DensityMeaInit,cellNumber,laneNumber, sample):
    mean = (DensityMeaInit[1]+DensityMeaInit[-2])/2
    std = 0.05*mean
    state = random.normal(mean, std, (sample,cellNumber))
    model = laneNumber * ones((sample, cellNumber))
    return state, model

####################################  end boundary conditions functions section   #######################################

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
    inflow = 1000
##  Model noise and measurement noise
    modelNoiseMean = 0.0
    modelNoiseStd = 5.0

    densityMeaMean = 0.0
    densityMeaStd = 13.5 #13.5

    speedMeaMean = -4.0
    speedMeaStd = 4.8 #4.4

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
    sample = 2500
    timeStep = 180
    pTran = 0.99
    pStay = 0.99
    pClear = 0.5

    PR = 4
    savefigure = True
    smoother = False
    lag = 2
    DeltaS = lag
    if smoother == False:
        lag=0

##  Discretization
    dx = Vmax*dt
    cellNumber = floor(length/dx)
    dx = length/cellNumber
    Lambda = dt/dx

## Load Measurements and true
    
##    directoryLoad = '/Users/Ren/Dropbox/SourceCode/Test/Corsim/SourceData/'
##    directorySave = '/Users/Ren/Dropbox/SourceCode/Test/Corsim/MMPF/Result_occ/'

    directoryLoad = os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir))+'/CORSIM filter factor/'
    directorySave = os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir))+'/test/MMPF/Result_inflow/'


    densityTrue = load(directoryLoad+'TrueState/'+'TrueDensityL3B1F'+str(inflow)+'S65R50.npy')
#    speedTrue = load(directoryLoad+'TrueSpeed.npy')
    densityMea = densityTrue
    speedMea = load(directoryLoad+'Measurements/'+'meaSpeed'+str(PR)+'R0L3B1F'+str(inflow)+'S65R50.npy')
#    parameterTrue = load(directoryLoad+'TrueParameter.npy')



## Creat array to save results
    estimatedState = zeros((timeStep, cellNumber))
    estimatedModel = zeros((timeStep, cellNumber))

## Initialization
    state, model = init_state(densityMea, cellNumber, laneNumber, sample)
    speed = Vmax*ones((sample,cellNumber))
    weight = init_weight(sample)

    estimatedState[0] = average(state, axis=0)
    estimatedModel[0] = laneNumber*ones(cellNumber)

##########################################################################################################################
## Algorithm Start
    start_time = time.time()

##    for counter in range(45,135):
    for counter in range(1,timeStep):
        print 'this is interation', counter
        
        if smoother == False or (timeStep-counter)<5:
            DeltaS = 1
        for i in range(DeltaS):

## measurements for the current time step
            densityMeaStep = zeros(cellNumber)
            densityMeaStep[1] = densityMea[counter+i, 1]
            densityMeaStep[-2] = densityMea[counter+i, -2]
            speedMeaStep = speedMea[counter+i]
## boundary conditions
            boundaryLeft = left_boundary_statistic(laneNumber, Vmax, inflow)
            boundaryRight = right_boundary_statistic(laneNumber, Vmax)

## regime transition
# no secondary incident case:
#            model = regime_transition_nb_normal(model, sample, pTran, laneNumber, cellNumber)
# secondary incident allowed:
            model = regime_transition_nb_second(model, sample, pTran, pStay, pClear, laneNumber, cellNumber)

# forward prediction
            for j in range(int(sample)):
                state[j] = ctm_lq(state[j],model[j], Lambda, boundaryLeft, boundaryRight, rhoc, rhom, Vmax, modelNoiseMean, modelNoiseStd, rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber)
                speed[j] = vel_lq(state[j], rhoc, rhom, Vmax, model[j], rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber)
                
            if i==0:
                stateCopy = state.copy()
                modelCopy = model.copy()

            likelihood = compute_particle_likelihood(state, speed, cellNumber, sample, densityMeaStep, speedMeaStep, densityMeaMean, densityMeaStd, speedMeaMean, speedMeaStd)
            weight = update_weight(likelihood, weight)
        state, model = resampling(stateCopy, modelCopy, sample, cellNumber, weight)
        weight = init_weight(sample)
        estimatedState[counter] = average(state, axis=0)
        estimatedModel[counter] = average(model, axis=0)
        print estimatedModel[counter]
        print 'run time for this iteration', time.time() - start_time
        


    save( directorySave + str(PR)+'lag'+str(lag)+'estDen'+str(inflow),estimatedState)
    save( directorySave + str(PR)+'lag'+str(lag)+'estModel'+str(inflow),estimatedModel)



    plot_density(densityTrue[0:180], savefigure, 'trueDensityPlot'+str(inflow), directorySave)
    plot_density(estimatedState, savefigure, str(PR)+'lag'+str(lag)+'estDen'+str(inflow), directorySave)
    plot_parameter(estimatedModel, savefigure, str(PR)+'lag'+str(lag)+'EstModel'+str(inflow), directorySave)
