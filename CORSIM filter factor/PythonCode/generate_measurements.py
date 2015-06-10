## This file generates the density and speed measurements 

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from numpy import *
from random import choice
from random import sample
from copy import deepcopy
import sys
import os

def build_veh_trajectory(keyData):
## construct the vehicle trajectory dictionary
## key: vehID
## value: [TimeIndex, VehPosition, Speed]
        vehTrajectory = dict()
        for i in range(len(keyData[:,0])):
                if str(int(keyData[i,1])) not in vehTrajectory:
                        vehTrajectory[str(int(keyData[i,1]))] = keyData[i,[0,2,3]]
                else:
                        vehTrajectory[str(int(keyData[i,1]))] = vstack((vehTrajectory[str(int(keyData[i,1]))], keyData[i,[0,2,3]]))
        return vehTrajectory

def build_data_time(keyData):
## construct the data time dictionary
## key: TimeIndex
## value: [VehID, VehPosition, Speed]
        dataTime = dict()
        for i in range(len(keyData[:,0])):
                if str(int(keyData[i,0])) not in dataTime:
                        dataTime[str(int(keyData[i,0]))] = keyData[i,[1,2,3]]
                else:
                        dataTime[str(int(keyData[i,0]))] = vstack((dataTime[str(int(keyData[i,0]))], keyData[i,[1,2,3]]))
        return dataTime

def collect_mea_density(vehTrajectory,dataTime, dt, dx, cellNumber, simulationTime):
## get the traffic density for all cells by collecting occupancy data
        timeStep = int(simulationTime/(dt))
        occupancyCount = zeros((timeStep, cellNumber))
        occupancy = zeros((timeStep, cellNumber))
        cellDensity = zeros((timeStep, cellNumber))
        for timeID in dataTime:
                timeIndex = int(int(timeID)/dt)
                try:
                        for j in range(len(dataTime[timeID][:,0])):
                                locationID = int(dataTime[timeID][j,1]/(dx*5280))
                                sensorCenterPoint = locationID*dx*5280 + dx*5280/2
                                if (sensorCenterPoint-10)<dataTime[timeID][j,1]<(sensorCenterPoint+10):
                                        occupancyCount[timeIndex, locationID] = occupancyCount[timeIndex, locationID]+1
                except IndexError:
                        locationID = int(dataTime[timeID][0,1]/(dx*5280))
                        sensorCenterPoint = locationID*dx*5280 + dx*5280/2
                        if (sensorCenterPoint-10)<dataTime[timeID][0,1]<(sensorCenterPoint+10):
                                occupancyCount[timeIndex, locationID] = occupancyCount[timeIndex, locationID]+1                                
        occupancy = occupancyCount/3/dt
        cellDensity = occupancy*5280.0/(14+6)*3
        return cellDensity


def collect_mea_speed(vehTrajectory, dataTime, dt, dx, cellNumber, simulationTime, penetrationRate):
## collect the speed measurements with the given penetration rate 
        timeStep = int(simulationTime/(dt))
        dictCellSpeed = dict()
        
        for timeID in dataTime:
                if mod(int(timeID),dt) == 0:
                        timeIndex = int(int(timeID)/dt)
                        try:
                                for j in range(len(dataTime[timeID][:,0])):
                                        locationID = int(dataTime[timeID][j,1]/(dx*5280))
                                        if locationID < 11.0: ## vehicle outside the domain
                                                if random.random()<(penetrationRate/100.0):
                                                        if (timeIndex, locationID) not in dictCellSpeed:
                                                                dictCellSpeed[(timeIndex, locationID)] = [dataTime[timeID][j,2]]
                                                        else:
                                                                dictCellSpeed[(timeIndex, locationID)].append(dataTime[timeID][j,2])
                        except IndexError:
                                if random.random()<(penetrationRate/100.0):
                                        locationID = int(dataTime[timeID][1]/(dx*5280))
                                        if (timeIndex, locationID) not in dictCellSpeed:
                                                dictCellSpeed[(timeIndex, locationID)] = [dataTime[timeID][2]]
                                        else:
                                                dictCellSpeed[(timeIndex, locationID)].append(dataTime[timeID][2])
        cellSpeed = zeros((timeStep, cellNumber))
        for element in dictCellSpeed:
                cellSpeed[element[0], element[1]] = average(dictCellSpeed[element])
        return cellSpeed
                                
def plot(cellDensity, cellSpeed, saveDirectory, fileName):
        plt.rc('xtick',labelsize=30)
        plt.rc('ytick',labelsize=30)
        plt.imshow(cellDensity,aspect='auto',origin='lower',interpolation='nearest')
        plt.ylabel('Time Step',fontsize=30)
        plt.clim(0.0, 500)
        plt.xlabel('Cell Number',fontsize=30)
        plt.colorbar()
        plt.savefig(saveDirectory+ 'meaDensity'+ str(penetrationRate) + fileName+'.pdf', bbox_inches='tight')
#        plt.show()
        plt.clf()

        cmap=plt.cm.jet_r
        plt.rc('xtick',labelsize=30)
        plt.rc('ytick',labelsize=30)
        plt.imshow(cellSpeed,aspect='auto',origin='lower',interpolation='nearest',cmap=cmap)
        plt.clim(0,80)
        plt.ylabel('Time Step',fontsize=30)
        plt.xlabel('Cell Number',fontsize=30)
        plt.colorbar()
        plt.savefig(saveDirectory+ 'meaSpeed'+ str(penetrationRate) + fileName+'.pdf', bbox_inches='tight')
#        plt.show()
        plt.clf()









if __name__=='__main__':



        loadDirectory = os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/ExtractedKeyData/'
        ## the directory where the output from this function will be saved at
        saveDirectory = os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Measurements/'  



        fileList = ['L3B1F5000S65R50', 'L3B1F4000S65R50', 'L3B1F3000S65R50', 'L3B1F2000S65R50', 'L3B1F1000S65R50']
##        fileList = ['L3B1F6000S65R50']
        fileNumber = len(fileList)

        ## discretization
        Vmax=65.0
        length=4.0
        dt = 20.0
        dx = Vmax*dt/3600.0
        cellNumber = int(length/dx)
        dx = length/cellNumber

        ## simulation parameter
        simulationTime = 3600
        penetrationRange = [4]
        randomSample = 1
        
        
        fileCount = 0
        for fileName in fileList:
                fileCount = fileCount+1
                print 'This is file', fileCount, 'out of', fileNumber

                keyData = load(loadDirectory+fileName+'key.npy')
                dataTime = build_data_time(keyData)
                vehTrajectory = build_veh_trajectory(keyData)
                cellDensityMea = collect_mea_density(vehTrajectory,dataTime, dt, dx, cellNumber, simulationTime)
                save(saveDirectory + 'meaDensity' + fileName, cellDensityMea)
                for penetrationRate in penetrationRange:
                        for rand in range(randomSample):                
                                cellSpeedMea = collect_mea_speed(vehTrajectory, dataTime, dt, dx, cellNumber, simulationTime, penetrationRate)
                                save(saveDirectory + 'meaSpeed'+ str(penetrationRate)+ 'R' +str(rand) + fileName, cellSpeedMea)
                        plot(cellDensityMea, cellSpeedMea, saveDirectory, fileName)




                
##
##
##
##
##
##
##
##
##
##
##                
##
##
##def build_veh_trajectory(keyData):
#### construct the vehicle trajectory dictionary
#### key: vehID
#### value: [TimeIndex, VehPosition, Speed]
##        vehTrajectory = dict()
##        for i in range(len(keyData[:,0])):
##                if str(int(keyData[i,1])) not in vehTrajectory:
##                        vehTrajectory[str(int(keyData[i,1]))] = keyData[i,[0,2,3]]
##                else:
##                        vehTrajectory[str(int(keyData[i,1]))] = vstack((vehTrajectory[str(int(keyData[i,1]))], keyData[i,[0,2,3]]))
##        return vehTrajectory
##
##def build_data_time(keyData):
#### construct the data time dictionary
#### key: TimeIndex
#### value: [VehID, VehPosition, Speed]
##        dataTime = dict()
##        for i in range(len(keyData[:,0])):
##                if str(int(keyData[i,0])) not in dataTime:
##                        dataTime[str(int(keyData[i,0]))] = keyData[i,[1,2,3]]
##                else:
##                        dataTime[str(int(keyData[i,0]))] = vstack((dataTime[str(int(keyData[i,0]))], keyData[i,[1,2,3]]))
##        return dataTime
##
##def collect_flow_speed_occupancy(vehTrajectory,dataTime, dt, cellNumber, simulationTime):
##        timeStep = int(simulationTime/(dt))
##
#### get the traffic density for all cells by collecting occupancy data
##        occupancyCount = zeros((timeStep, cellNumber))
##        occupancy = zeros((timeStep, cellNumber))
##        cellDensity = zeros((timeStep, cellNumber))
##        for timeID in dataTime:
##                timeIndex = int(int(timeID)/dt)
##                try:
##                        for j in range(len(dataTime[timeID][:,0])):
##                                locationID = int(dataTime[timeID][j,1]/(dx*5280))
##                                sensorCenterPoint = locationID*dx*5280 + dx*5280/2
##                                if (sensorCenterPoint-10)<dataTime[timeID][j,1]<(sensorCenterPoint+10):
##                                        occupancyCount[timeIndex, locationID] = occupancyCount[timeIndex, locationID]+1
##                except IndexError:
##                        locationID = int(dataTime[timeID][0,1]/(dx*5280))
##                        sensorCenterPoint = locationID*dx*5280 + dx*5280/2
##                        if (sensorCenterPoint-10)<dataTime[timeID][0,1]<(sensorCenterPoint+10):
##                                occupancyCount[timeIndex, locationID] = occupancyCount[timeIndex, locationID]+1                                
##        occupancy = occupancyCount/3/dt
##        cellDensity = occupancy*5280.0/(14+6)*3
##
#### get the speed anf flow data from the trajectory data
##        dictFlowSpeed = dict() ## key: intervalIndex value: list of speeds
##        for vehID in vehTrajectory:
##                try:
##                        vehAppearNumber = len(vehTrajectory[vehID][:,1])
##                        for i in range(int(vehAppearNumber-1)):
##                                timeIndex = int(vehTrajectory[vehID][i,0]/dt)
##                                locationIndex = int(vehTrajectory[vehID][i,1]/(dx*5280))
##                                cellCenterPoint = locationIndex*dx*5280 + dx*5280/2
##                                speedCenter = vehTrajectory[vehID][i,2]
##
##                                if vehTrajectory[vehID][i,1] <= cellCenterPoint and vehTrajectory[vehID][i+1,1]>cellCenterPoint:
##                                        if (timeIndex, locationIndex) not in dictFlowSpeed:
##                                                dictFlowSpeed[(timeIndex, locationIndex)] = [speedCenter]
##                                        else:
##                                                dictFlowSpeed[(timeIndex, locationIndex)].append(speedCenter)
##                except IndexError:
##                        pass
##        cellSpeed = zeros((timeStep, cellNumber))
##        for element in dictFlowSpeed:
##                cellSpeed[element[0], element[1]] = average(dictFlowSpeed[element])
##
##
##        plt.rc('xtick',labelsize=30)
##        plt.rc('ytick',labelsize=30)
##        plt.imshow(cellDensity,aspect='auto',origin='lower',interpolation='nearest')
##        plt.ylabel('Time Step',fontsize=30)
##        plt.clim(0.0, 500)
##        plt.xlabel('Cell Number',fontsize=30)
##        plt.colorbar()
##        plt.savefig('TrueDensity.pdf', bbox_inches='tight')
##        plt.show()
##        plt.clf()
##
##
##        cmap=plt.cm.jet_r
##        plt.rc('xtick',labelsize=30)
##        plt.rc('ytick',labelsize=30)
##        plt.imshow(cellSpeed,aspect='auto',origin='lower',interpolation='nearest',cmap=cmap)
##        plt.clim(0,80)
##        plt.ylabel('Time Step',fontsize=30)
##        plt.xlabel('Cell Number',fontsize=30)
##        plt.colorbar()
##        plt.savefig('TrueSpeed.pdf', bbox_inches='tight')
##        plt.show()
##        plt.clf()
##        return cellDensity, cellSpeed
##
##
##def collect_flow_speed_occupancy_penetration(vehTrajectory,dataTime, dt, cellNumber, simulationTime, penetrationRate):
##        timeStep = int(simulationTime/(dt))
##
#### get the traffic density for all cells by collecting occupancy data
##        occupancyCount = zeros((timeStep, cellNumber))
##        occupancy = zeros((timeStep, cellNumber))
##        cellDensity = zeros((timeStep, cellNumber))
##        for timeID in dataTime:
##                timeIndex = int(int(timeID)/dt)
##                try:
##                        for j in range(len(dataTime[timeID][:,0])):
##                                locationID = int(dataTime[timeID][j,1]/(dx*5280))
##                                sensorCenterPoint = locationID*dx*5280 + dx*5280/2
##                                if (sensorCenterPoint-10)<dataTime[timeID][j,1]<(sensorCenterPoint+10):
##                                        occupancyCount[timeIndex, locationID] = occupancyCount[timeIndex, locationID]+1
##                except IndexError:
##                        locationID = int(dataTime[timeID][0,1]/(dx*5280))
##                        sensorCenterPoint = locationID*dx*5280 + dx*5280/2
##                        if (sensorCenterPoint-10)<dataTime[timeID][0,1]<(sensorCenterPoint+10):
##                                occupancyCount[timeIndex, locationID] = occupancyCount[timeIndex, locationID]+1                                
##        occupancy = occupancyCount/3/dt
##        cellDensity = occupancy*5280.0/(14+6)*3
##
#### get the speed anf flow data from the trajectory data
##        dictFlowSpeed = dict() ## key: intervalIndex value: list of speeds
##        for vehID in vehTrajectory:
##                try:
##                        vehAppearNumber = len(vehTrajectory[vehID][:,1])
##                        for i in range(int(vehAppearNumber-1)):
##                                timeIndex = int(vehTrajectory[vehID][i,0]/dt)
##                                locationIndex = int(vehTrajectory[vehID][i,1]/(dx*5280))
##                                cellCenterPoint = locationIndex*dx*5280 + dx*5280/2
##                                speedCenter = vehTrajectory[vehID][i,2]
##
##                                if vehTrajectory[vehID][i,1] <= cellCenterPoint and vehTrajectory[vehID][i+1,1]>cellCenterPoint:
##                                        if random.random()<(penetrationRate/100.0):
##                                                if (timeIndex, locationIndex) not in dictFlowSpeed:
##                                                        dictFlowSpeed[(timeIndex, locationIndex)] = [speedCenter]
##                                                else:
##                                                        dictFlowSpeed[(timeIndex, locationIndex)].append(speedCenter)
##                except IndexError:
##                        pass
##        cellSpeed = zeros((timeStep, cellNumber))
##        for element in dictFlowSpeed:
##                cellSpeed[element[0], element[1]] = average(dictFlowSpeed[element])
##
##
##        plt.rc('xtick',labelsize=30)
##        plt.rc('ytick',labelsize=30)
##        plt.imshow(cellDensity[0:180],aspect='auto',origin='lower',interpolation='nearest')
##        plt.ylabel('Time Step',fontsize=30)
##        plt.clim(0.0, 500)
##        plt.xlabel('Cell Number',fontsize=30)
##        plt.colorbar()
##        plt.savefig('meaDensity'+str(penetrationRate)+'.pdf', bbox_inches='tight')
##        plt.show()
##        plt.clf()
##
##
##        cmap=plt.cm.jet_r
##        plt.rc('xtick',labelsize=30)
##        plt.rc('ytick',labelsize=30)
##        plt.imshow(cellSpeed[0:180],aspect='auto',origin='lower',interpolation='nearest',cmap=cmap)
##        plt.clim(0,80)
##        plt.ylabel('Time Step',fontsize=30)
##        plt.xlabel('Cell Number',fontsize=30)
##        plt.colorbar()
##        plt.savefig('meaSpeed'+str(penetrationRate)+'.pdf', bbox_inches='tight')
##        plt.show()
##        plt.clf()
##        return cellDensity, cellSpeed
##
##
##
##if __name__=='__main__':
##
##        filename='3lane1blockS65F7000'
##        keyData=load(filename+'key.npy')
##
#### discretization
##        Vmax=65.0
##        length=4.0
##        dt = 20.0
##        dx = Vmax*dt/3600.0
##        cellNumber = int(length/dx)
##        dx = length/cellNumber
##
#### simulation parameters 
##        simulationTime=3800
##        penetrationRate = 10
##
##        dataTime = build_data_time(keyData)
##        vehTrajectory = build_veh_trajectory(keyData)
###        cellDensityTrue, cellSpeedTrue = collect_flow_speed_occupancy(vehTrajectory,dataTime, dt, cellNumber, simulationTime)
##        cellDensity, cellSpeed = collect_flow_speed_occupancy_penetration(vehTrajectory,dataTime, dt, cellNumber, simulationTime, penetrationRate)
##
##        save('trueDensity', cellDensity[0:180])
##        save('meaSpeed'+str(penetrationRate), cellSpeed[0:180])
##        
##
##
##
##
##        
####                                        
##
##
##
##
##
##                                
##                                if 990<dataTime[timeID][j,1]<1010:
##                                        collectedOccupancy[timeIndex] = collectedOccupancy[timeIndex]+1
##
##                except IndexError:
##                        pass
##
##
##
##        
##
##        
##        dataFlowSpeed = dict() ## key: intervalIndex value: list of speeds
##        collectedData = zeros((timeInterval-1,3)) # [density speed flow]
##        for vehID in vehTrajectory:
##                try:
##                        vehAppearNumber = len(vehTrajectory[vehID][:,1])
##                        if vehAppearNumber > 1:
##                                for i in range(int(vehAppearNumber-1)):
##                                        if vehTrajectory[vehID][i,1] <= 1000 and vehTrajectory[vehID][i+1,1]>1000:
##                                                timeSensor =  vehTrajectory[vehID][i,0]
##                                                speedSensor = vehTrajectory[vehID][i,2]
##                                                if timeSensor > deltaT:
##                                                        intervalIndex = int(timeSensor/deltaT)-1
##                                                        if str(intervalIndex) not in dataFlowSpeed:
##                                                                dataFlowSpeed[str(intervalIndex)] = [speedSensor]
##                                                        else:
##                                                                dataFlowSpeed[str(intervalIndex)].append(speedSensor)
##                except IndexError:
##                        pass
##
##        collectedOccupancy = zeros(timeInterval-1)        
##        for timeID in dataTime:
##                if int(timeID) >= deltaT:
##                        timeIndex = int(int(timeID)/deltaT)-1
##                        try:
##                                for j in range(len(dataTime[timeID][:,0])):
##                                        if 990<dataTime[timeID][j,1]<1010:
##                                                collectedOccupancy[timeIndex] = collectedOccupancy[timeIndex]+1
##
##                        except IndexError:
##                                pass
##        collectedOccupancy = collectedOccupancy/3/deltaT
##        for k in dataFlowSpeed:
##                flow = len(dataFlowSpeed[k])*(simulationTime/deltaT)
##                speed = average(dataFlowSpeed[k])
##                collectedData[k,0]=collectedOccupancy[int(k)]*5280/(14+6)*3
##                collectedData[k,1]=speed
##                collectedData[k,2]=flow
##        return collectedData
##
##














####filename='3lane1blockS65F7000'
####keyData=load(filename+'key.npy')
####
####
###### Define functions
####
###### Construct the DataTime dictionary. Key: time. Values: VehID/Position/Speed
####def buildDataTime(keyData):
####        DataTime=dict()
####        for i in range(int(len(keyData[:,0]))):
####                if str(int(keyData[i,0])) not in DataTime:
####                        DataTime[str(int(keyData[i,0]))]=keyData[i,[1,2,3]]
####                else:
####                        DataTime[str(int(keyData[i,0]))]=vstack((DataTime[str(int(keyData[i,0]))],keyData[i,[1,2,3]]))
####        return DataTime
####
###### Construct the VehTrajectory dictionaray. Key: VehID. Values: Time/Position/Speed
####def buildVehTrajectory(keyData):
####        VehTrajectory=dict()
####        for i in range(int(len(keyData[:,1]))):
####                if str(int(keyData[i,1])) not in VehTrajectory:
####                        VehTrajectory[str(int(keyData[i,1]))]=keyData[i,[0,2,3]]
####                else:
####                        VehTrajectory[str(int(keyData[i,1]))]=vstack((VehTrajectory[str(int(keyData[i,1]))],keyData[i,[0,2,3]]))
####        return VehTrajectory
####
###### Identify GPS vehicles
###### Construct DataTimeSelected dictionary, DataTimeSelected is a subset of DataTime, when GPS vehicles present
####def getGPSDataTime(headway, DataTime):       
####        DataTimeSelected=dict()
####        GPSlist=[]
####        for i in list(DataTime):
####                if mod(int(i),headway)==0:
####                        try:
####                                j=DataTime[i][:,1].argmin()
####                                GPSlist.append(DataTime[i][j,0])
####                        except IndexError:
####                                GPSlist.append(DataTime[i][0])
####        for i in list(DataTime):
####                try:
####                        for j in range(len(DataTime[i][:,0])):
####                                if DataTime[i][j,0] in GPSlist:
####                                        if i not in DataTimeSelected:
####                                                DataTimeSelected[i]=DataTime[i][j,:]
####                                        else:
####                                                DataTimeSelected[i]=vstack((DataTimeSelected[i],DataTime[i][j,:]))
####                except IndexError:
####                        if DataTime[i][0] in GPSlist:
####                                DataTimeSelected[i]=DataTime[i][:]
####        return GPSlist, DataTimeSelected
####
####
###### Calculate cellDensity
####def calculateCellDensity(SimulationTime, TimeInterval, cellNumber, keyData, spaceStep,length):       
####        NoVeh=zeros((int(SimulationTime/TimeInterval),cellNumber))
####        for k in range(int(len(keyData[:,0]))):
####                if mod(int(keyData[k,0]),TimeInterval)==0:
####                        if keyData[k,2]>=length*5280:
####                                pass
####                        else:
####                                j=int(keyData[k,2]/(spaceStep*5280))
####                                NoVeh[int(keyData[k,0])/TimeInterval,j]=NoVeh[int(keyData[k,0])/TimeInterval,j]+1
####        cellDensity=NoVeh/spaceStep
####        return NoVeh, cellDensity
####
####
###### Construct TimeCell dictionary Key-Value(Key)-Value time->cells->[vehID,speed]
####def calculateTimeCell(DataTime,spaceStep,length):
####        TimeCell=dict()
####        for i in list(DataTime):
####                CellVehSpeed=dict()
####                try:
####                        for j in range(len(DataTime[i][:,1])):
####                                position=DataTime[i][j,1]
####                                if position>=length*5280:
####                                        pass
####                                else:                                       
####                                        cell=int(position/(spaceStep*5280))
####                                        if str(cell) not in CellVehSpeed:
####                                                CellVehSpeed[str(cell)]=DataTime[i][j,[0,2]]
####                                        else:
####                                                CellVehSpeed[str(cell)]=vstack((CellVehSpeed[str(cell)],DataTime[i][j,[0,2]]))
####                except IndexError:
####                        position=DataTime[i][1]
####                        if position>=length*5280:
####                                pass
####                        else:
####                                cell=int(position/(spaceStep*5280))
####                                if str(cell) not in CellVehSpeed:
####                                        CellVehSpeed[str(cell)]=DataTime[i][[0,2]]
####                                else:
####                                        CellVehSpeed[str(cell)]=vstack((CellVehSpeed[str(cell)],DataTime[i][[0,2]]))                       
####                TimeCell[i]=CellVehSpeed
####        return TimeCell
####
###### Calculate cellSpeed
####def calculateCellSpeed(SimulationTime, TimeInterval, cellNumber, TimeCell):       
####        cellSpeed=zeros((int(SimulationTime/TimeInterval),cellNumber))
####        for i in list(TimeCell):
####                if mod(int(i),TimeInterval)==0:
####                        for j in list(TimeCell[i]):
####                                try:
####                                        cellSpeed[int(i)/TimeInterval,int(j)]=average(TimeCell[i][j][:,1])
####                                except IndexError:
####                                        cellSpeed[int(i)/TimeInterval,int(j)]=TimeCell[i][j][1]
####        return cellSpeed
####
####
####
####
###### Here starts the main code
####
####
####Vmax=65.0
####length=4.0
####
####TimeInterval=20
####timeStep=TimeInterval/3600.0
####SimulationTime=3800
####
####spaceStep=Vmax*timeStep
####cellNumber=floor(length/spaceStep)
####spaceStep=length/cellNumber ## this is unitLength in terms of mile
####
###### This is the penetration rate value
####headway=20
####
####
####DataTime=buildDataTime(keyData)
####print 'Finish DataTime'
####
####VehTrajectory=buildVehTrajectory(keyData)
####print 'Finish vehTrajectory'
####
####GPSlist, GPSDataTime=getGPSDataTime(headway, DataTime)
####print 'Finish GPSDataTime'
####
####NoVeh, cellDensity=calculateCellDensity(SimulationTime, TimeInterval, cellNumber, keyData, spaceStep,length)
####print 'Finish cellDensity'
####
####TimeCell=calculateTimeCell(GPSDataTime, spaceStep, length)
####print 'Finish TimeCell'
####
####
####cellSpeed=calculateCellSpeed(SimulationTime, TimeInterval, cellNumber, TimeCell)
####print 'Finish cellSpeed'
####
####
####
####save('TrueDensity', cellDensity)
####save('RawSpeed'+str(headway), cellSpeed)
####
####
######## Here ends the main code
####
####vehID=list(VehTrajectory)
####
####plt.hold(True)
####for i in vehID:
####        try:
####                plt.plot(VehTrajectory[i][:,0],VehTrajectory[i][:,1])
####        except IndexError:
####                ## vehicle trajectory may contain single dot, for example, a vehicle enters at the end of the simulation
####                pass
####plt.xlabel('Time')
####plt.ylabel('Space')
####plt.savefig('Raw'+filename+str(headway)+'diagram.pdf')
####plt.hold(False)
####plt.clf()
####
####
####imgplot1=plt.imshow(cellDensity,aspect='auto',origin='lower',interpolation='none')
####plt.ylabel('Time Step')
####plt.xlabel('Cell Number')
####plt.colorbar()
######plt.show()
####plt.savefig('TrueDensity.pdf')
####plt.clf()
####
####
####imgplot2=plt.imshow(cellSpeed,aspect='auto',origin='lower',interpolation='none')
####plt.ylabel('Time Step')
####plt.xlabel('Cell Number')
####plt.colorbar()
######plt.show()
####plt.savefig('Raw'+filename+str(headway)+'Speed.pdf')
####plt.clf()
