import matplotlib.pyplot as plt
from numpy import *
from random import choice
from copy import deepcopy




def build_system_models_nb(cellNumber, laneNumber, maxLaneBlocked):
##  This function constructs the model set of the system assuming no incidents occur at the two boundary cells
##  Used for CORSIM implementation
    modelNumber = (cellNumber-4)*maxLaneBlocked+1
    modelSet = laneNumber*ones((modelNumber,cellNumber))
    i=1
    for j in range(2,int(cellNumber-2)):
        for k in range(1,maxLaneBlocked+1):
            modelSet[i,j]=laneNumber-k
            i=i+1
    return modelSet, modelNumber

def build_system_models(cellNumber, laneNumber, maxLaneBlocked):
##  This function constructs the model set of the system
    modelNumber = cellNumber*maxLaneBlocked+1
    modelSet = laneNumber*ones((modelNumber,cellNumber))
    i=1
    for j in range(0,int(cellNumber)):
        for k in range(1,maxLaneBlocked+1):
            modelSet[i,j]=laneNumber-k
            i=i+1
    return modelSet, modelNumber

def init_model_probability_nb(previousMode, pTran, modelNumber, laneNumber, maxLaneBlocked):
##  Initialize the model probability for all models in the system, given the previous selected model, assuming no incidents at the two boundary cells
##  Used for CORSIM implementation
    modelProbability = ones(modelNumber)*(1-pTran)/(modelNumber-1)
    modelProbability[0] = pTran
    for i in range(2,len(previousMode)-2):
        for k in range(1,maxLaneBlocked+1):
            if previousMode[i] == (laneNumber-k):
                modelProbability[0]=(1-pTran)/(modelNumber-1)
                modelProbability[k+maxLaneBlocked*(i-2)]=pTran
    return modelProbability

def init_model_probability_nb_back(previousMode, pTran, modelNumber, laneNumber, maxLaneBlocked):
##  Initialize the model probability for all models in the system, given the previous selected model, assuming no incidents at the two boundary cells
##  Used for CORSIM implementation, a higher probability to transition back to no incident.
    modelProbability = ones(modelNumber)*(1-pTran)/(modelNumber-1)
    modelProbability[0] = pTran
    for i in range(2,len(previousMode)-2):
        for k in range(1,maxLaneBlocked+1):
            if previousMode[i] == (laneNumber-k):
                modelProbability[0]=(1-(pTran-0.1))/(modelNumber-1)
                modelProbability[k+maxLaneBlocked*(i-2)]=pTran-0.1
    return modelProbability

def init_model_probability(previousModel, pTran, modelNumber, laneNumber, maxLaneBlocked):
##  Initialize the model probability for all models in the system, given the previous selected model
    modelProbability = ones(modelNumber)*(1-pTran)/(modelNumber-1)
    modelProbability[0] = pTran
    for i in range(0,len(previousModel)):
        for k in range(1,maxLaneBlocked+1):
            if previousModel[i] == (laneNumber-k):
                modelProbability[0]=(1-pTran)/(modelNumber-1)
                modelProbability[k+maxLaneBlocked*i]=pTran
    return modelProbability


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
    

def compute_KalmanGain(cellDensity, cellSpeed, cellNumber, sample, densityMea, speedMea, densityMeaStd, speedMeaStd):
## compute the kalman gain
    H=matrix(build_operatorH(cellNumber,densityMea, speedMea))
    meaNumber = sum(H)
    densityMean = matrix(average(cellDensity,axis=0)).T
    densitySpeed = hstack((cellDensity, cellSpeed))
    densitySpeedMean = matrix(average(densitySpeed,axis=0)).T
# create covariance matrix
    covMatrixPH = zeros((cellNumber, meaNumber))
    covMatrixHPH = zeros((meaNumber, meaNumber))
    covMatrixR = zeros((meaNumber, meaNumber))
    for i in range(int(sample)):
        covMatrixPH = covMatrixPH + (matrix(cellDensity[i]).T-densityMean)* ((H*(matrix(densitySpeed[i]).T)-H*densitySpeedMean).T)
        covMatrixHPH = covMatrixHPH +(H*(matrix(densitySpeed[i]).T)-H*densitySpeedMean)*((H*(matrix(densitySpeed[i]).T)-H*densitySpeedMean).T)
    covMatrixPH = covMatrixPH/(sample-1)
    covMatrixHPH = covMatrixHPH/(sample-1)
    for j in range(int(meaNumber)):
        if j<2:
            covMatrixR[j,j] = densityMeaStd*densityMeaStd
        else:
            covMatrixR[j,j] = speedMeaStd*speedMeaStd
    kalmanGain = covMatrixPH*((covMatrixHPH+covMatrixR).I)
    return kalmanGain
        

def kalman_state_update(cellDensity, cellSpeed, cellNumber, sample, densityMea, speedMea, kalmanGain, densityMeaMean, densityMeaStd, speedMeaMean, speedMeaStd):
# update the system state
    H=matrix(build_operatorH(cellNumber,densityMea, speedMea))
    meaNumber = sum(H)
    densitySpeed = hstack((cellDensity, cellSpeed))
    allMea = hstack((densityMea, speedMea))
    meaNoise=zeros(meaNumber)
    for i in range(0,int(meaNumber)):
        if i<sum(densityMea !=0):
            meaNoise[i] = random.normal(densityMeaMean, densityMeaStd)
        else:
            meaNoise[i] = random.normal(speedMeaMean, speedMeaStd)
    for j in range(int(sample)):
        cellDensity[j] = cellDensity[j] + array(matrix(kalmanGain*(H*(matrix(allMea).T)-H*(matrix(densitySpeed[j]).T)+matrix(meaNoise).T)).T)
    return cellDensity            
    


def compute_model_likelihood(cellDensity, cellSpeed, cellNumber, densityMea, speedMea, densityMeaMean, densityMeaStd, speedMeaMean, speedMeaStd):
# compute the likelihood for a model
    H = build_operatorH(cellNumber, densityMea, speedMea)
    allMea = hstack((densityMea, speedMea))
    estimatedDensitySpeedMean = hstack((average(cellDensity,axis=0),average(cellSpeed,axis=0)))
    diff = (matrix(H)*matrix(allMea).T-matrix(H)*matrix(estimatedDensitySpeedMean).T)
    modelLikelihood= 1.0 
    for i in range(0,int(sum(H))):
        if i<sum(densityMea !=0):
            modelLikelihood=modelLikelihood*1/(densityMeaStd*sqrt(2*pi))*exp(-(diff[i]-densityMeaMean)*(diff[i]-densityMeaMean)/(2*densityMeaStd*densityMeaStd))
        else:
            modelLikelihood=modelLikelihood*1/(speedMeaStd*sqrt(2*pi))*exp(-(diff[i]-speedMeaMean)*(diff[i]-speedMeaMean)/(2*speedMeaStd*speedMeaStd))
    return modelLikelihood





























    
    
    















