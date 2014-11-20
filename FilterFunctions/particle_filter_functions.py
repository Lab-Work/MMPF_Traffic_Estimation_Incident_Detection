import matplotlib.pyplot as plt
from numpy import *
from random import choice
from copy import deepcopy

def init_weight(sample):
## particle weight update function
    weight = 1.0/sample*ones(sample)
    return weight

def regime_transition_nb_normal(model, sample, pTran, laneNumber, cellNumber):
## Regime transition, assuming no incident at the two boundary cells, one incident at a time.
    for i in range(sample):
            if min(model[i]) == laneNumber:
                    if random.random() > pTran:
                            incidentPosition=choice(range(2,int(cellNumber-2)))
                            model[i][incidentPosition]=choice(range(0,laneNumber))
            else:
                    if random.random()>pTran:
                            model[i]=laneNumber*ones(cellNumber)
    return model

def regime_transition_normal(model, sample, pTran, laneNumber, cellNumber):
## Regime transition, one incident at a time.
    for i in range(sample):
            if min(model[i]) == laneNumber:
                    if random.random() > pTran:
                            incidentPosition=choice(range(int(cellNumber)))
                            model[i][incidentPosition]=choice(range(0,laneNumber))
            else:
                    if random.random()>pTran:
                            model[i]=laneNumber*ones(cellNumber)
    return model

def regime_transition_nb_second(model, sample, pTran, pStay, pClear, laneNumber, cellNumber):
## Regime transition, assuming no incident at the two boundary cells, a second incident may occur on the upsteam of the existing incident 
    for i in range(sample):
            if min(model[i]) == laneNumber:
                    if random.random() > pTran:
                            incidentPosition=choice(range(2,int(cellNumber-2)))
                            model[i][incidentPosition]=choice(range(0,laneNumber))
            else:
                incidentNumber=0
                for j in range(len(model[i])):
                    if model[i][j] != laneNumber:
                        incidentNumber = incidentNumber+1
                        if incidentNumber == 1:
                            incidentLocationOne = j
                            incidentSeverityOne = model[i][j]
                        elif incidentNumber == 2:
                            incidentLocationTwo = j
                            incidentSeverityTwo = model[i][j]
                        else:
                            'something is wrong with the incident transition model, please check'
                if incidentNumber == 1:
                    if random.random() > pStay:
                        if random.random() > pClear:
                            try:
                                incidentPositionSecond = choice(range(2,int(incidentLocationOne)))
                                incidentSeverityList = range(int(laneNumber))
                                indexRemove = incidentSeverityList.index(incidentSeverityOne)
                                del incidentSeverityList[indexRemove]
                                incidentSeveritySecond = choice(incidentSeverityList)
                            except IndexError:
                                incidentPositionSecond = 2
                                incidentSeverityList = range(int(laneNumber))
                                indexRemove = incidentSeverityList.index(incidentSeverityOne)
                                del incidentSeverityList[indexRemove]
                                incidentSeveritySecond = choice(incidentSeverityList)
                            model[i][incidentPositionSecond] = incidentSeveritySecond
                        else:
                            model[i] = laneNumber*ones(cellNumber)
                elif incidentNumber == 2:
                    if random.random() >pStay:
                            incidentLocationClear = choice([incidentLocationOne, incidentLocationTwo])
                            model[i][incidentLocationClear] = laneNumber
    return model



def build_operatorH(cellNumber, densityMea, speedMea):
##  Construct the linear operator for the nonlinear observation equation
##  H is m by 2n, m is the number of measurements, n is the number of cells
    mea=hstack((densityMea,speedMea))
    meaExist= mea!=0
    meaNumber=int(sum(meaExist))
    H=zeros((meaNumber,int(2*cellNumber)))
    if meaNumber==0:
        print 'there is no measurements at current time step'
    else:
        k=0
        for i in range(int(2*cellNumber)):
            if mea[i] != 0:
                H[k,i]=1
                k=k+1
    return H

def compute_particle_likelihood(cellDensity, cellSpeed, cellNumber, sample, densityMea, speedMea, densityMeaMean, densityMeaStd, speedMeaMean, speedMeaStd):
## Calculate the likelihood of each particle
    likelihood = zeros(int(sample))
    H=matrix(build_operatorH(cellNumber,densityMea, speedMea))
    for j in range(sample):
        allMea = hstack((densityMea, speedMea))
        estimatedDensitySpeed = hstack((cellDensity[j], cellSpeed[j]))
        diff = (matrix(H)*matrix(allMea).T-matrix(H)*matrix(estimatedDensitySpeed).T)
        modelLikelihood= 1.0 
        for i in range(0,int(sum(H))):
            if i<sum(densityMea !=0):
                modelLikelihood=modelLikelihood*1.0/(densityMeaStd*sqrt(2*pi))*exp(-(diff[i]-densityMeaMean)*(diff[i]-densityMeaMean)/(2*densityMeaStd*densityMeaStd))
            else:
                modelLikelihood=modelLikelihood*1.0/(speedMeaStd*sqrt(2*pi))*exp(-(diff[i]-speedMeaMean)*(diff[i]-speedMeaMean)/(2*speedMeaStd*speedMeaStd))
        likelihood[j] = modelLikelihood
    likelihood = likelihood/sum(likelihood)
    return likelihood

def update_weight(likelihood, weight):
## particle weight update function 
    weight=likelihood*weight
    weight=weight/sum(weight)
    return weight

def resampling(state, model, sample, cellNumber, weight):
    Cum=cumsum(weight)
    stateCopy = state.copy()
    modelCopy = model.copy()
    
    step=1.0/sample
    i=1
    u1=random.uniform(0,step)
    for j in range(sample):
        uj = u1+step*(j-1)
        while uj>Cum[i]:
            i=i+1
        state[j]=stateCopy[i]
        model[j]=modelCopy[i]
    return (state, model)



## Description of parameters

# cellNumber : the number of cells
# densityMea : density measurements
# speedMea   : speed measurements
# EstDensity : estimated density
# EstSpeed   : estimated speed
# sample     : the number of samples for the particle filter algorithm
# meanDensity: the mean of the noise model of density measurements
# meanSpeed  : the mean of the noise model of speed measurements
# stdDensity : the std of the noise model of density measurements
# stdSpeed   : the std of the noise model of speed measurements

# EnNumber   : the number of ensembles for the Kalman filter algorithm
# DensityMea : the density measurements (dimension n)
# SpeedMea   : the speed measurements (dimension n)

# Ren Wang
# UIUC, June 2014

##
##
##def convert_measurement(cellNumber, densityMea, speedMea):
#### Convert the 2n dimensional measurement vector to the m dimension vector
#### m is the number of measurements available, n is the number of cells
##    MeaDS=hstack((densityMea,speedMea))
##    H=build_operatorH(cellNumber, densityMea, speedMea)
##    Mea=matrix(MeaDS)*matrix(H).T
##    return array(Mea)
##
##def convert_state(EstDensity, EstSpeed, cellNumber, densityMea, speedMea):
#### Nonlinear observation operator which matches the 2n dimensional states to m dimensional measurements
##    H=build_operatorH(cellNumber, densityMea, speedMea)
##    stateMea=hstack((EstDensity, EstSpeed))
##    stateMea=matrix(stateMea)*matrix(H).T
##    return array(stateMea)
##
##def compute_normalize_particle_likelihood(EstDensity, EstSpeed, densityMea, speedMea, sample, cellNumber, meanDensity, meanSpeed, stdDensity, stdSpeed):
#### Calculate the likelihood of each particle 
##    MeaDSexist = densityMea != 0
##    MeaSPexist = speedMea != 0
##    MeaDSNumber = int(sum(MeaDSexist))
##    MeaSPNumber = int(sum(MeaSPexist))
### get the location of density and speed measurements
##    MeaDindex = where(MeaDSexist==1)
##    MeaSindex = where(MeaSPexist==1)
##
##    DiffDensity = EstDensity-densityMea
##    DiffSpeed = EstSpeed-speedMea
##    likelihoodDensity=1/(stdDensity*sqrt(2*pi))*exp(-(DiffDensity-meanDensity)*(DiffDensity-meanDensity)/(2*stdDensity*stdDensity))
##    likelihoodSpeed=1/(stdSpeed*sqrt(2*pi))*exp(-(DiffSpeed-meanSpeed)*(DiffSpeed-meanSpeed)/(2*stdSpeed*stdSpeed))
##
##    likelihoodD=ones((sample, cellNumber))
##    likelihoodS=ones((sample, cellNumber))
##    likelihoodD[:,MeaDindex]=likelihoodDensity[:,MeaDindex]
##    likelihoodS[:,MeaSindex]=likelihoodSpeed[:,MeaSindex]
##    
##    likelihood=likelihoodD*likelihoodS
####    likelihood=likelihood.chip(min=1e-20)
##    likelihood_particle=prod(likelihood,axis=1)
##    Sum=sum(likelihood_particle)
####  print Sum
##    Likelihood_particle_normalize=likelihood_particle/Sum
##    return Likelihood_particle_normalize
##
##def compute_particle_likelihood(EstDensity, EstSpeed, densityMea, speedMea, sample, cellNumber, meanDensity, meanSpeed, stdDensity, stdSpeed):
#### Calculate the likelihood of each particle 
##    MeaDSexist = densityMea != 0
##    MeaSPexist = speedMea != 0
##    MeaDSNumber = int(sum(MeaDSexist))
##    MeaSPNumber = int(sum(MeaSPexist))
### get the location of density and speed measurements
##    MeaDindex = where(MeaDSexist==1)
##    MeaSindex = where(MeaSPexist==1)
##
##    DiffDensity = EstDensity-densityMea
##    DiffSpeed = EstSpeed-speedMea
##    likelihoodDensity=1/(stdDensity*sqrt(2*pi))*exp(-(DiffDensity-meanDensity)*(DiffDensity-meanDensity)/(2*stdDensity*stdDensity))
##    likelihoodSpeed=1/(stdSpeed*sqrt(2*pi))*exp(-(DiffSpeed-meanSpeed)*(DiffSpeed-meanSpeed)/(2*stdSpeed*stdSpeed))
##
##    likelihoodD=ones((sample, cellNumber))
##    likelihoodS=ones((sample, cellNumber))
##    likelihoodD[:,MeaDindex]=likelihoodDensity[:,MeaDindex]
##    likelihoodS[:,MeaSindex]=likelihoodSpeed[:,MeaSindex]
##    
##    likelihood=likelihoodD*likelihoodS
####    likelihood=likelihood.chip(min=1e-20)
##    likelihood_particle=prod(likelihood,axis=1)
##    return likelihood_particle
##
##def update_weight(likelihood, Weight):
#### particle weight update function 
##    Weight=likelihood*Weight
##    Weight=Weight/sum(Weight)
##    return Weight
##
##def update_weight_notNormalized(likelihood, Weight):
#### particle weight update function 
##    Weight=likelihood*Weight
##    return Weight
##
##
##def resampling(EstDensity, Model, sample, cellNumber, Weight):
##    Weight=Weight/sum(Weight)
##    Cum=cumsum(Weight)
##    
##    EstDensityStore = deepcopy(EstDensity)
##    ModelStore = deepcopy(Model)
##
##    
##    step=1.0/sample
##    i=1
##    u1=random.uniform(0,step)
##    for j in range(sample):
##        uj = u1+step*(j-1)
##        while uj>Cum[i]:
##            i=i+1
##        EstDensity[j]=EstDensityStore[i]
##        Model[j]=ModelStore[i]
##    return (EstDensity, Model)
##
##def regime_transition(Model, sample, pTran, laneNumber, cellNumber):
#### Regime transition 
##    for i in range(sample):
##            if min(Model[i]) == laneNumber:
##                    if random.random() > pTran:
##                            incidentPosition=int(cellNumber*random.random())
##                            Model[i][incidentPosition]=choice(range(0,laneNumber))
##            else:
##                    if random.random()>pTran:
##                            Model[i]=laneNumber*ones(cellNumber)
##    return Model
##
##def build_system_mode(cellNumber, laneNumber, maxLaneBlocked):
#### build all possible modes of the system
##    ModeNumber = cellNumber*maxLaneBlocked+1
##    Mode = laneNumber*ones((ModeNumber,cellNumber))
##    i=1
##    for j in range(int(cellNumber)):
##        for k in range(1,maxLaneBlocked+1):
##            Mode[i,j]=laneNumber-k
##            i=i+1
##    return Mode, ModeNumber
##
##def build_system_mode_nb(cellNumber, laneNumber, maxLaneBlocked):
#### build all possible modes of the system, no incident on boundary cells
##    ModeNumber = (cellNumber-4)*maxLaneBlocked+1
##    Mode = laneNumber*ones((ModeNumber,cellNumber))
##    i=1
##    for j in range(2,int(cellNumber-2)):
##        for k in range(1,maxLaneBlocked+1):
##            Mode[i,j]=laneNumber-k
##            i=i+1
##    return Mode, ModeNumber
##
##def init_model_probability(previousMode, pTran, ModelNumber):
##    modelProbability = ones(ModelNumber)*(1-pTran)/(ModelNumber-1)
##    modelProbability[0] = pTran
##    for i in range(2,len(previousMode)-2):
##        if previousMode[i]==2:
##            modelProbability[0]=(1-pTran)/(ModelNumber-1)
##            modelProbability[1+3*(i-2)]=pTran
##        elif previousMode[i]==1:
##            modelProbability[0]=(1-pTran)/(ModelNumber-1)
##            modelProbability[2+3*(i-2)]=pTran                
##        elif previousMode[i]==0:
##            modelProbability[0]=(1-pTran)/(ModelNumber-1)
##            modelProbability[3+3*(i-2)]=pTran
##    return modelProbability
##
##def build_operatorH(cellNumber, densityMea, speedMea):
#### The linear operator for the nonlinear observation equation
### H is m by 2n, m is the number of measurements, n is the number of cells
##    MeaDS=hstack((densityMea,speedMea))
##    MeaExist= MeaDS!=0
##    meaNumber=int(sum(MeaExist))
##    H=zeros((meaNumber,int(2*cellNumber)))
##    if meaNumber==0:
##        print 'there is no measurements at current time step'
##    else:
##        k=0
##        for i in range(int(2*cellNumber)):
##            if MeaDS[i] != 0:
##                H[k,i]=1
##                k=k+1
##    return H
##    
##def compute_Cov(cellDensity, cellNumber, EnNumber):
### compute the covariance matrix
##    CovMatrix=zeros((cellNumber, cellNumber))
##    for i in range(0,EnNumber):
##        CovMatrix=CovMatrix+matrix(cellDensity[i,:]-average(cellDensity,axis=0)).T*matrix(cellDensity[i,:]-average(cellDensity,axis=0))
##    CovMatrix=CovMatrix/(EnNumber-1)
##    return CovMatrix
##
##def compute_KalmanGain(cellDensity, cellSpeed, cellNumber,EnNumber, densityMea, speedMea, stdDensity, stdSpeed):
### compute the kalman gain
##    H=build_operatorH(cellNumber,densityMea, speedMea)
##    StatePredict=hstack((cellDensity, cellSpeed))
##    StatePredictMean=average(StatePredict,axis=0)
##    CovMatrixHPH=zeros((sum(H), sum(H)))
##    CovMatrixPH=zeros((cellNumber,sum(H)))
##    CovMatrixR=zeros((sum(H), sum(H)))
##    for i in range(0,EnNumber):
##        CovMatrixHPH=CovMatrixHPH+matrix(matrix(H)*matrix(StatePredict[i,:]).T-matrix(H)*matrix(StatePredictMean).T).T*matrix(matrix(H)*matrix(StatePredict[i,:]).T-matrix(H)*matrix(StatePredictMean).T)
##        CovMatrixPH=CovMatrixPH+matrix(cellDensity[i,:]-average(cellDensity,axis=0)).T*matrix(matrix(H)*matrix(StatePredict[i,:]).T-matrix(H)*matrix(StatePredictMean).T).T
##    CovMatrixHPH=CovMatrixHPH/(EnNumber-1)
##    CovMatrixPH=CovMatrixPH/(EnNumber-1)
##    for j in range(0,int(sum(H))):
##        if j<2:
##            CovMatrixR[j,j]=stdDensity*stdDensity
##        else:
##            CovMatrixR[j,j]=stdSpeed*stdSpeed
##    kalmanGain=matrix(CovMatrixPH)*matrix(CovMatrixHPH+CovMatrixR).I
##    return kalmanGain
##
##
##def state_estimation(cellDensity,cellSpeed,cellNumber,EnNumber, densityMea, speedMea, kalmanGain, stdDensity, stdSpeed):
### estimate the state
##    StatePredict=hstack((cellDensity, cellSpeed))
##    H=build_operatorH(cellNumber, densityMea, speedMea)
##    TrueMea=hstack((densityMea, speedMea))
##    meaNoise=zeros(sum(H))
##    for i in range(0,int(sum(H))):
##        if i<sum(densityMea !=0):
##            meaNoise[i]=stdDensity*stdDensity
##        else:
##            meaNoise[i]=stdSpeed*stdSpeed
##    for i in range(0,EnNumber):
##        cellDensity[i,:]=cellDensity[i,:]+array((kalmanGain*(matrix(H)*matrix(TrueMea).T-matrix(H)*matrix(StatePredict[i,:]).T+matrix(meaNoise).T)).T)[0]
##    return cellDensity
##
##def compute_mode_likelihood(cellDensity, cellSpeed, cellNumber, densityMea, speedMea, meanDensity, meanSpeed, stdDensity, stdSpeed):
### compute the likelihood for a mode
##    H=build_operatorH(cellNumber, densityMea, speedMea)
##    TrueMea=hstack((densityMea, speedMea))
##    EstimationStateMean=hstack((average(cellDensity,axis=0),average(cellSpeed,axis=0)))
##    Diff=(matrix(H)*matrix(TrueMea).T-matrix(H)*matrix(EstimationStateMean).T)
##    modeLikelihood= 1.0 
##    for i in range(0,int(sum(H))):
##        if i<sum(densityMea !=0):
##            modeLikelihood=modeLikelihood*1/(stdDensity*sqrt(2*pi))*exp(-(Diff[i]-meanDensity)*(Diff[i]-meanDensity)/(2*stdDensity*stdDensity))
##        else:
##            modeLikelihood=modeLikelihood*1/(stdSpeed*sqrt(2*pi))*exp(-(Diff[i]-meanSpeed)*(Diff[i]-meanSpeed)/(2*stdSpeed*stdSpeed))
##    return modeLikelihood
##
##def imm_factor(previousMode, pTran, ModeNumber, laneNumber, maxLaneBlocked):
### build the weight for each mode based on the transition probability. 
##    adjustFactor = ones(ModeNumber)*(1-pTran)/(ModeNumber-1)
##    adjustFactor[0]=pTran
##    for i in range(len(previousMode)):
##        for j in range(1,int(maxLaneBlocked+1)):
##            if previousMode[i] == laneNumber-j:
##                adjustFactor[0] = (1-pTran)/(ModeNumber-1)
##                index=j+maxLaneBlocked*i
##                adjustFactor[index]=pTran
##    return adjustFactor
##
##
##
##
##
##
##
##
##










































    
    
    















