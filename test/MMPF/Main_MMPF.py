import matplotlib.pyplot as plt
from numpy import *
from random import choice
from copy import deepcopy
import sys
##sys.path.append('/Users/Ren/Dropbox/SourceCode/FilterFunctions/')
##sys.path.append('/Users/Ren/Dropbox/SourceCode/TrafficModel/')
sys.path.append('/Users/Ren/Dropbox/SourceCode/FilterFunctions/')
sys.path.append('/Users/Ren/Dropbox/SourceCode/TrafficModel/')
from traffic_models_lq import *
from particle_filter_functions import *
import time

####################################  statictical boundary conditions   #######################################

def left_boundary_statistic(laneNumber, Vmax):
    lbDensity = (5900.0+ random.normal(0,150))/Vmax
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
##  Model noise and measurement noise
    modelNoiseMean = 5.0
    modelNoiseStd = 20.0

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
    sample = 250
    timeStep = 180
    pTran = 0.99
    pStay = 0.99
    pClear = 0.5

    PR = 4
    savefigure = True
    smoother = False
    lag = 4
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

    directoryLoad = '/Users/Ren/Dropbox/SourceCode/CORSIM filter factor/'
    directorySave = '/Users/Ren/Dropbox/SourceCode/test/MMPF/Result_occ/'


    densityTrue = load(directoryLoad+'TrueState/'+'TrueDensityL3B1F6000S65R50.npy')
#    speedTrue = load(directoryLoad+'TrueSpeed.npy')
    densityMea = densityTrue
    speedMea = load(directoryLoad+'Measurements/'+'meaSpeed'+str(PR)+'R0L3B1F6000S65R50.npy')
#    parameterTrue = load(directoryLoad+'TrueParameter.npy')



## Creat array to save results
    estimatedState = zeros((timeStep, cellNumber))
    estimatedModel = zeros((timeStep, cellNumber))
    estimatedStateAll = zeros((timeStep, sample, cellNumber))
    estimatedModelAll = zeros((timeStep, sample, cellNumber))

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
            boundaryLeft = left_boundary_statistic(laneNumber, Vmax)
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
        estimatedStateAll[counter] = state
        estimatedModelAll[counter] = model
        print estimatedModel[counter]
        print 'run time for this iteration', time.time() - start_time
        

#    stateError = estimatedState - densityTrue[0:180]
#    modelError = estimatedModel - parameterTrue[0:180]

    print "run time is", time.time() - start_time
#    print 'density error is', average(abs(stateError))
#    print 'model error is', sum(abs(modelError))

    save( directorySave + str(PR)+'lag'+str(lag)+'estDen',estimatedState)
    save( directorySave + str(PR)+'lag'+str(lag)+'estModel',estimatedModel)
    save( directorySave + str(PR)+'lag'+str(lag)+'estDenAll',estimatedStateAll)
    save( directorySave + str(PR)+'lag'+str(lag)+'estModelAll',estimatedModelAll)
    
#    save( directorySave + str(headway)+'lag'+str(lag)+'stateError',stateError)
#    save( directorySave + str(headway)+'lag'+str(lag)+'modelError',modelError)


    plot_density(densityTrue[0:180], savefigure, 'trueDensityPlot', directorySave)
    plot_density(estimatedState, savefigure, str(PR)+'lag'+str(lag)+'estDen', directorySave)
    plot_parameter(estimatedModel, savefigure, str(PR)+'lag'+str(lag)+'EstModel', directorySave)
#    plot_parameter(parameterTrue[0:180], savefigure, 'trueModelPlot', directorySave)

##
##    
##
##plot_density(DensityTrue[:,1:-1], savefigure, 'TrueDensityPlot'+'lag'+str(lag), directory)
##plot_speed(SpeedTrue[:,1:-1], savefigure, 'TrueSpeedPlot', directory)
##
##plot_density(EstimatedDensity1st, savefigure, str(headway)+'lag'+str(lag)+'EstDen1st', directory)
##plot_density(EstimatedDensity2nd, savefigure, str(headway)+'lag'+str(lag)+'EstDen2nd', directory)
##plot_speed(EstimatedSpeed1st, savefigure, str(headway)+'lag'+str(lag)+'EstSpeed1st', directory)
##plot_speed(EstimatedSpeed2nd, savefigure, str(headway)+'lag'+str(lag)+'EstSpeed2nd', directory)
##plot_parameter(EstimatedModel1st, savefigure, str(headway)+'lag'+str(lag)+'EstModel1st', directory)
##plot_parameter(EstimatedModel2nd, savefigure, str(headway)+'lag'+str(lag)+'EstModel2nd', directory)
##plot_density(DensityError1st, savefigure, str(headway)+'lag'+str(lag)+'ErrDen1st', directory)
##plot_density(DensityError2nd, savefigure, str(headway)+'lag'+str(lag)+'ErrDen2nd', directory)
##plot_speed(SpeedError1st, savefigure, str(headway)+'lag'+str(lag)+'ErrSpeed1st', directory)
##plot_speed(SpeedError2nd, savefigure, str(headway)+'lag'+str(lag)+'ErrSpeed2nd', directory)
##
##
##



#### Creat array to save results
##    EstimatedDensity1st = zeros((timeStep,cellNumber))
##    EstimatedDensity2nd = zeros((timeStep,cellNumber))
##    EstimatedSpeed1st = zeros((timeStep,cellNumber))
##    EstimatedSpeed2nd = zeros((timeStep,cellNumber))
##    EstimatedModel1st = zeros((timeStep,cellNumber))
##    EstimatedModel2nd = zeros((timeStep,cellNumber))
##
#### Initialization
##
### 1st model
##    state_1st, model_1st = init_state_1st(DensityMea,cellNumber,laneNumber, sample)
##    speed_1st = init_speed(sample, cellNumber)
##    weight_1st = init_weight(sample)
##    EstimatedDensity1st[0] = average(state_1st, axis=0)
##    EstimatedSpeed1st[0] = average(speed_1st, axis=0)
##    EstimatedModel1st[0] = laneNumber*ones(cellNumber)
##
### 2nd model
##    state_2nd, model_2nd = init_state_2nd(DensityMea,cellNumber,laneNumber, sample) 
##    speed_2nd = init_speed(sample, cellNumber)
##    weight_2nd = init_weight(sample)
##    EstimatedDensity2nd[0] = average(state_2nd[:,0,:], axis=0)
##    EstimatedSpeed2nd[0] = average(speed_2nd, axis=0)
##    EstimatedModel2nd[0] = laneNumber*ones(cellNumber)
##    
##
#### Algorithm Start
##    start_time = time.time()
##    for counter in range(1,timeStep):
##        print 'this is interation', counter
##                
##        if Smoother == False or (timeStep-counter)<5:
##            lag = 1
##        for i in range(lag):
##
#### Attach the boundary conditions
##            lbc_2nd = left_boundary(laneNumber,DensityTrue, SpeedTrue, counter+i)
##            rbc_2nd = right_boundary(laneNumber,DensityTrue, SpeedTrue, counter+i)
##            lbc_1st = lbc_2nd[0],lbc_2nd[2]
##            rbc_1st = rbc_2nd[0],rbc_2nd[2]
##
##            wl,wr = boundary_qq(lbc_2nd,rbc_2nd,rhoc1,rhoc2,rhom1,rhom2,Vmax,rhom_tuide, NewtonIteration)
##            
##            lbc_2nd = lbc_2nd[0], wl, lbc_2nd[2] 
##            rbc_2nd = rbc_2nd[0], wr, rbc_2nd[2]
##
##
#### regime transition
##            model_1st = regime_transition(model_1st, sample, ptran, laneNumber, cellNumber)
##            model_2nd = regime_transition(model_2nd, sample, ptran, laneNumber, cellNumber)
##
##            for j in range(int(sample)):
##
##                state_1st[j] = ctm_qq(state_1st[j],model_1st[j], Lambda, lbc_1st, rbc_1st, rhoc, rhom, Vmax, rhom_tuide,ModelNoiseMean, ModelNoiseStd)
##                speed_1st[j] = vel_qq(state_1st[j], model_1st[j]*rhoc,model_1st[j]*rhom,Vmax,model_1st[j]*rhom_tuide)
##
##                state_2nd[j] = ctm_qq_2nd(state_2nd[j], model_2nd[j],Lambda, lbc_2nd, rbc_2nd, rhoc1, rhoc2, rhom1, rhom2, Vmax, rhom_tuide,ModelNoiseMean, ModelNoiseStd)
##                speed_2nd[j] = vel_qq(state_2nd[j,0], model_2nd[j]*(rhoc1*rhoc2/(state_2nd[j,1]*rhoc2+(1-state_2nd[j,1])*rhoc1)), model_2nd[j]*(rhom1*rhom2/(state_2nd[j,1]*rhom2+(1-state_2nd[j,1])*rhom1)), Vmax, model_2nd[j]*rhom_tuide)
##
##            if i == 0:
##
##                state_1st_copy = state_1st.copy()
##                model_1st_copy = model_1st.copy()
##
##                state_2nd_copy = state_2nd.copy()
##                model_2nd_copy = model_2nd.copy()
##
##            likelihood_1st = compute_normalize_particle_likelihood(state_1st, speed_1st, zeros(cellNumber), SpeedMea[counter+i], sample, cellNumber, DensityMeaMean, SpeedMeaMean, DensityMeaStd, SpeedMeaStd)
##            likelihood_2nd = compute_normalize_particle_likelihood(state_2nd[:,0,:], speed_2nd, zeros(cellNumber), SpeedMea[counter+i], sample, cellNumber, DensityMeaMean, SpeedMeaMean, DensityMeaStd, SpeedMeaStd)
##
##            weight_1st = update_weight(likelihood_1st, weight_1st)
##            weight_2nd = update_weight(likelihood_2nd, weight_2nd)
##
##        state_1st, model_1st = resampling(state_1st_copy,model_1st_copy, sample, cellNumber, weight_1st)
##        state_2nd, model_2nd = resampling(state_2nd_copy, model_2nd_copy, sample, cellNumber, weight_2nd)
##
##        weight_1st = init_weight(sample)
##        weight_2nd = init_weight(sample)
##
##        EstimatedDensity1st[counter] = average(state_1st, axis=0)
##        EstimatedDensity2nd[counter] = average(state_2nd[:,0,:], axis=0)
##        EstimatedSpeed1st[counter] = average(speed_1st, axis=0)
##        EstimatedSpeed2nd[counter] = average(speed_2nd, axis=0)
##        EstimatedModel1st[counter] = average(model_1st, axis=0)
##        EstimatedModel2nd[counter] = average(model_2nd, axis=0)
##
##        print 'model prob for 1st model is:', average(model_1st, axis=0)
##        print 'model prob for 2nd model is:', average(model_2nd, axis=0)
##
##    DensityError1st = EstimatedDensity1st-DensityTrue[:,1:-1]
##    DensityError2nd = EstimatedDensity2nd-DensityTrue[:,1:-1]
##    SpeedError1st = EstimatedSpeed1st-SpeedTrue[:,1:-1]
##    SpeedError2nd = EstimatedSpeed2nd-SpeedTrue[:,1:-1]
##
##    print "run time is", time.time() - start_time
##    print 'density error is', average(abs(DensityError2nd))
##
##    save( directory +'MMPF/' +'Result/' + str(headway)+'lag'+str(lag)+'EstDen1st',EstimatedDensity1st)
##    save( directory +'MMPF/' +'Result/' + str(headway)+'lag'+str(lag)+'EstDen2nd',EstimatedDensity2nd)
##    save( directory +'MMPF/' +'Result/' + str(headway)+'lag'+str(lag)+'EstSpeed1st',EstimatedSpeed1st)
##    save( directory +'MMPF/' +'Result/' + str(headway)+'lag'+str(lag)+'EstSpeed2nd',EstimatedSpeed2nd)
##    save( directory +'MMPF/' +'Result/' + str(headway)+'lag'+str(lag)+'EstModel1st',EstimatedModel1st)
##    save( directory +'MMPF/' +'Result/' + str(headway)+'lag'+str(lag)+'EstModel2nd',EstimatedModel2nd)
##
##
##
##plot_density(DensityTrue[:,1:-1], savefigure, 'TrueDensityPlot', directory)
##plot_speed(SpeedTrue[:,1:-1], savefigure, 'TrueSpeedPlot', directory)
##
##plot_density(EstimatedDensity1st, savefigure, str(headway)+'lag'+str(lag)+'EstDen1st', directory)
##plot_density(EstimatedDensity2nd, savefigure, str(headway)+'lag'+str(lag)+'EstDen2nd', directory)
##plot_speed(EstimatedSpeed1st, savefigure, str(headway)+'lag'+str(lag)+'EstSpeed1st', directory)
##plot_speed(EstimatedSpeed2nd, savefigure, str(headway)+'lag'+str(lag)+'EstSpeed2nd', directory)
##plot_parameter(EstimatedModel1st, savefigure, str(headway)+'lag'+str(lag)+'EstModel1st', directory)
##plot_parameter(EstimatedModel2nd, savefigure, str(headway)+'lag'+str(lag)+'EstModel2nd', directory)
##plot_density(DensityError1st, savefigure, str(headway)+'lag'+str(lag)+'ErrDen1st', directory)
##plot_density(DensityError2nd, savefigure, str(headway)+'lag'+str(lag)+'ErrDen2nd', directory)
##plot_speed(SpeedError1st, savefigure, str(headway)+'lag'+str(lag)+'ErrSpeed1st', directory)
##plot_speed(SpeedError2nd, savefigure, str(headway)+'lag'+str(lag)+'ErrSpeed2nd', directory)
##


