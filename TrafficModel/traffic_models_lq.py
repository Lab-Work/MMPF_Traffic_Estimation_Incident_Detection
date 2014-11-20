import matplotlib.pyplot as plt
from numpy import *
from random import choice
from copy import deepcopy




########################################################################################################################################################
########################################################################################################################################################

def ctm_lq(rhoVector, modelVector, Lambda, bdl, bdr, rhoc, rhom, Vmax, ModelNoiseMean, ModelNoiseStd, rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber, inflow = -10, outflow = -10):
        
## The linear quadratic(lq) model is a CTM with a linear quadratic fundamental diagram
        
## Description of parameters
#  rhoVector : system state at time n
#  modelVector: system model at time n, a vector which specifies the number of lanes open at each cell
#  Lambda =dt/dx    : discretization parameter in the Godunov scheme
#  bdl=[rhol,modell] : upstream density and model
#  bdr=[rhor,modelr] : downstream density and model
#  rhoc: critical density
#  rhom: maximum density
#  Vmax: maximum speed
#  ModelNoiseMean: mean value of the CTM noise
#  ModelNoiseStd: std value of the CTM noise
#  rhocOneIncident: critical density for incident cell, when one lane is blocked
#  rhocTwoIncident: critical density for incident cell, when two lanes are blocked
#  VmaxIncident: maximum speed for incident cell
#  laneNumber: total number of lanes of the road
#  inflow : inflow to the link (used for the network solver, not used for the MMPF paper)
#  outflow: outflow from the link (used for the network solver, not used for the MMPF paper)

# Ren Wang


## the boundary condition
        rhoLeftGhost = bdl[0]
        rhoRightGhost = bdr[0]
        modelLeftGhost = bdl[1]
        modelRightGhost = bdr[1]

## forward prediction
        rhoVectorUpdated = rhoVector.copy()
        for i in range(len(rhoVector)):
                if i==0:
                        if inflow >= 0:
                                rhoVectorUpdated[i] = rhoVector[i] + Lambda*(inflow-\
                                                                 min(sending_lq(rhoVector[i], rhoc, rhom, Vmax, modelVector[i], rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber), receiving_lq(rhoVector[i+1], rhoc, rhom, Vmax, modelVector[i+1], rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber)))
                        else:                
                                rhoVectorUpdated[i] = rhoVector[i] + Lambda*(min(sending_lq(rhoLeftGhost, rhoc, rhom, Vmax, modelLeftGhost, rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber), receiving_lq(rhoVector[i], rhoc, rhom, Vmax, modelVector[i], rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber))-\
                                                                 min(sending_lq(rhoVector[i], rhoc, rhom, Vmax, modelVector[i], rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber), receiving_lq(rhoVector[i+1], rhoc, rhom, Vmax, modelVector[i+1], rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber)))
                elif i<int(len(rhoVector)-1):
                        rhoVectorUpdated[i] = rhoVector[i] + Lambda*(min(sending_lq(rhoVector[i-1], rhoc, rhom, Vmax, modelVector[i-1], rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber), receiving_lq(rhoVector[i], rhoc, rhom, Vmax, modelVector[i], rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber))-\
                                                         min(sending_lq(rhoVector[i], rhoc, rhom, Vmax, modelVector[i], rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber), receiving_lq(rhoVector[i+1], rhoc, rhom, Vmax, modelVector[i+1], rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber)))
                else:
                        if outflow >= 0:
                                rhoVectorUpdated[i] = rhoVector[i] + Lambda*(min(sending_lq(rhoVector[i-1], rhoc, rhom, Vmax, modelVector[i-1], rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber), receiving_lq(rhoVector[i], rhoc, rhom, Vmax, modelVector[i], rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber))-outflow)
                        else:
                                rhoVectorUpdated[i] = rhoVector[i] + Lambda*(min(sending_lq(rhoVector[i-1], rhoc, rhom, Vmax, modelVector[i-1], rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber), receiving_lq(rhoVector[i], rhoc, rhom, Vmax, modelVector[i], rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber))-\
                                                                 min(sending_lq(rhoVector[i], rhoc, rhom, Vmax, modelVector[i], rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber), receiving_lq(rhoRightGhost, rhoc, rhom, Vmax, modelRightGhost, rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber)))
        return rhoVectorUpdated+random.normal(ModelNoiseMean, ModelNoiseStd, (len(rhoVector)))

def sending_lq(rho, rhoc, rhom, Vmax, model, rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber):
        if model == laneNumber or model == 0.0:
                rhoc = model*rhoc
        elif model == 2.0:
                Vmax = VmaxIncident
                rhoc = model*rhocOneIncident
        elif model == 1.0:
                Vmax = VmaxIncident
                rhoc = model*rhocTwoIncident
        else:
                print 'there is a model not defined in the estimation model'
                
        if rhoc == 0 or rho<0:
                qSend = 0.0
        elif rho<rhoc:
                qSend = rho*Vmax
        else:
                qSend = rhoc*Vmax
        return qSend

def receiving_lq(rho, rhoc, rhom, Vmax, model, rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber):
        if model == laneNumber or model == 0.0:
                rhoc = model*rhoc
                rhom = model*rhom
        elif model == 2.0:
                Vmax = VmaxIncident
                rhoc = model*rhocOneIncident
                rhom = model*rhom
        elif model == 1.0:
                Vmax = VmaxIncident
                rhoc = model*rhocTwoIncident
                rhom = model*rhom
        else:
                print 'there is a model not defined in the estimation model'

        if rhoc ==0 or rho>rhom:
                qReceive = 0.0
        elif rho<rhoc:
                qReceive = rhoc*Vmax
        else:
                qmax = rhoc*Vmax
                a = qmax/(2*rhoc*rhom-rhom*rhom-rhoc*rhoc)
                b = -2*rhoc*qmax/(2*rhoc*rhom-rhom*rhom-rhoc*rhoc)
                c = qmax*(2*rhoc*rhom-rhom*rhom)/(2*rhoc*rhom-rhom*rhom-rhoc*rhoc)
                qReceive = a*rho*rho+b*rho+c
        return qReceive

def vel_lq(rhoVector, rhoc, rhom, Vmax, modelVector, rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber):
        VmaxCopy = Vmax
        speed = zeros(len(rhoVector))
        for i in range(len(rhoVector)):
                if modelVector[i] == laneNumber or modelVector[i] == 0.0:
                        Vmax = VmaxCopy
                        rhocCell = rhoc*modelVector[i]
                        rhomCell = rhom*modelVector[i]
                elif modelVector[i] == 2.0:
                        Vmax = VmaxIncident
                        rhocCell = rhocOneIncident*modelVector[i]
                        rhomCell = rhom*modelVector[i]
                elif modelVector[i] == 1.0:
                        Vmax = VmaxIncident
                        rhocCell = rhocTwoIncident*modelVector[i]
                        rhomCell = rhom*modelVector[i]
                        
                rhoCell = rhoVector[i]
                if rhocCell == 0 or rhoCell>rhomCell:
                        v = 0.0
                elif rhoCell<0:
                        v = Vmax
                elif rhoCell<rhocCell:
                        v = Vmax
                else:
                        qmax = rhocCell*Vmax
                        a = qmax/(2*rhocCell*rhomCell-rhomCell*rhomCell-rhocCell*rhocCell)
                        b = -2*rhocCell*qmax/(2*rhocCell*rhomCell-rhomCell*rhomCell-rhocCell*rhocCell)
                        c = qmax*(2*rhocCell*rhomCell-rhomCell*rhomCell)/(2*rhocCell*rhomCell-rhomCell*rhomCell-rhocCell*rhocCell)
                        v = (a*rhoCell*rhoCell+b*rhoCell+c)/rhoCell
                speed[i] = v
        return speed








###################################################### Test first order model ##########################################################################


##def left_boundary_statistic(laneNumber, Vmax):
##    lbDensity = (5900.0+ random.normal(0,250))/Vmax
##    lbModel = laneNumber
##    return lbDensity, lbModel
##
##def right_boundary_statistic(laneNumber, Vmax):
##    rbDensity = (5900.0+ random.normal(0,250))/Vmax
##    rbModel = laneNumber
##    return rbDensity, rbModel
##
##def plot_density(data):
##    plt.rc('xtick',labelsize=30)
##    plt.rc('ytick',labelsize=30)
##    plt.imshow(data,aspect='auto',origin='lower',interpolation='nearest')
##    plt.ylabel('Time Step',fontsize=30)
###    plt.clim(0.0, 560)
##    plt.xlabel('Cell Number',fontsize=30)
##    plt.colorbar()
##    plt.show()
##    plt.clf()
##
##
####  Model noise and measurement noise
##modelNoiseMean = 0.0
##modelNoiseStd = 1.0
##
##densityMeaMean = 0.0
##densityMeaStd = 7.5
##
##speedMeaMean = 0.0
##speedMeaStd = 4.3
##
####  Parameters for traffic models
##rhoc = 35.0
##rhom = 239.0    
##Vmax = 65.0
##rhom_tuide = 1303.0
##rhocOneIncident = 90.4# 115 #143.0
##rhocTwoIncident = 62.7
##VmaxIncident = 18.0#12.8#13.8
##rhom_tuideIncident = 30000.0
##
####  Parameters for simulation
##dt = 20.0 / 3600
##length = 4.0
##laneNumber = 3
##maxLaneBlocked = 3
##sample = 10
##timeStep = 180
##pTran = 0.99
##headway = 20
##savefigure = True
##
####  Discretization
##dx = Vmax*dt
##cellNumber = floor(length/dx)
##dx = length/cellNumber
##Lambda = dt/dx
##
####  state
##estimatedState = zeros((timeStep, cellNumber))
##estimatedSpeed = zeros((timeStep, cellNumber))
##state = zeros(cellNumber)
##model = laneNumber*ones(cellNumber)
##model[4]=3
##bdl=[90,3]
##bdr=[90,3]
##
##for i in range(int(timeStep)):
##        if i==60:
##                model[4]=2
##        if i==120:
##                model[4]=3
##        state=ctm_lq(state, model, Lambda, bdl, bdr, rhoc, rhom, Vmax, modelNoiseMean, modelNoiseStd, rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber)
##        estimatedSpeed[i]=vel_lq(state, rhoc, rhom, Vmax, model, rhocOneIncident, rhocTwoIncident, VmaxIncident, laneNumber)
##        estimatedState[i]=state.copy()
##        
##
##plot_density(estimatedState)
##plt.rc('xtick',labelsize=30)
##plt.rc('ytick',labelsize=30)
##plt.imshow(estimatedSpeed,aspect='auto',origin='lower',interpolation='nearest')
##plt.ylabel('Time Step',fontsize=30)
##plt.clim(0.0, 80)
##plt.xlabel('Cell Number',fontsize=30)
##plt.colorbar()
##plt.show()
##plt.clf()


###################################################### Test second order model ##########################################################################

##def plot_density(data):
##    plt.rc('xtick',labelsize=30)
##    plt.rc('ytick',labelsize=30)
##    plt.imshow(data,aspect='auto',origin='lower',interpolation='nearest')
##    plt.ylabel('Time Step',fontsize=30)
###    plt.clim(0.0, 560)
##    plt.xlabel('Cell Number',fontsize=30)
##    plt.colorbar()
##    plt.show()
##    plt.clf()
##
#### 2nd test
##dt=20.0/3600
##dx=64.8*dt
##bdl=[105,0.5,3]
##bdr=[105,0.5,3]
##Vmax=65
##
####state = zeros(cellNumber)
####model = laneNumber*ones(cellNumber)
##
##Lambda=dt/dx
##length = 4.0
##cellNumber = floor(length/dx)
##
##rhoVector = 10*ones(cellNumber)
##wVector = 0.5*ones(cellNumber)
##modelVector = 3.0* ones(cellNumber)
##
######  Parameters for traffic models
##rhoc1 = 37.0
##rhoc2 = 33.0
##rhom1 = 243.0
##rhom2 = 235.0
##Vmax = 65.0
##rhom_tuide = 1303.0
##rhoc1OneIncident = 92.4
##rhoc2OneIncident = 88.4
##rhoc1TwoIncident = 64.7
##rhoc2TwoIncident = 60.7
##VmaxIncident = 18.0
##rhom_tuideIncident = 30000.0
##laneNumber = 3.0
##ModelNoiseMean=0
##ModelNoiseStd=0.1
##timeStep=180
##
####  state
##estimatedState = zeros((timeStep, cellNumber))
##estimatedSpeed = zeros((timeStep, cellNumber))
##state = zeros(cellNumber)
##model = laneNumber*ones(cellNumber)
##model[4]=3
##
##for i in range(int(timeStep)):
##        if i==60:
##                modelVector[4]=2
##        if i==120:
##                modelVector[4]=3
##        rhoVector, wVector = ctm_2qq_incident(rhoVector, wVector, modelVector, Lambda, bdl, bdr,\
##                     rhoc1, rhoc2, rhom1, rhom2, Vmax, rhom_tuide,\
##                     rhoc1OneIncident, rhoc2OneIncident, rhoc1TwoIncident, rhoc2TwoIncident,\
##                     VmaxIncident, rhom_tuideIncident, laneNumber, ModelNoiseMean, ModelNoiseStd)
##        estimatedSpeed[i] = vel_2qq_incident(rhoVector, wVector, rhoc1, rhoc2, rhom1, rhom2, Vmax, modelVector, rhom_tuide, \
##                     rhoc1OneIncident, rhoc2OneIncident, rhoc1TwoIncident, rhoc2TwoIncident,\
##                     VmaxIncident, rhom_tuideIncident, laneNumber)
##
##
##        
##        estimatedState[i]=rhoVector.copy()
##        
##
##plot_density(estimatedState)
##plt.rc('xtick',labelsize=30)
##plt.rc('ytick',labelsize=30)
##plt.imshow(estimatedSpeed,aspect='auto',origin='lower',interpolation='nearest')
##plt.ylabel('Time Step',fontsize=30)
##plt.clim(0.0, 80)
##plt.xlabel('Cell Number',fontsize=30)
##plt.colorbar()
##plt.show()
##plt.clf()



#########################################################################################################################################################

