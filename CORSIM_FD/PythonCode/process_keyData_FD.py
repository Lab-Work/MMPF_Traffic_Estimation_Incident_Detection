## This function process the keyData file for fundamental diagram calibration
## keyData: [TimeIndex, VehID, VehPosition, Speed]

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from numpy import *
from random import choice
from random import sample
from copy import deepcopy


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

def collect_flow_speed_occupancy(vehTrajectory, dataTime, deltaT, simulationTime, sensorLocation, sensorLength):
## collect flow, speed and occupancy data from a sensor located at sensorLocation with length: sensorLength 

        timeInterval = simulationTime/deltaT
        dataFlowSpeed = dict() ## key: intervalIndex value: list of speeds
        collectedData = zeros((timeInterval-1,3)) # [density, speed, flow]
        for vehID in vehTrajectory:
                try:
                        vehAppearNumber = len(vehTrajectory[vehID][:,1])
                        if vehAppearNumber > 1:
                                for i in range(int(vehAppearNumber-1)):
                                        if vehTrajectory[vehID][i,1] <= sensorLocation and vehTrajectory[vehID][i+1,1]>sensorLocation:
                                                timeSensor =  vehTrajectory[vehID][i,0]
                                                speedSensor = vehTrajectory[vehID][i,2]
                                                if timeSensor > deltaT:
                                                        intervalIndex = int(timeSensor/deltaT)-1
                                                        if str(intervalIndex) not in dataFlowSpeed:
                                                                dataFlowSpeed[str(intervalIndex)] = [speedSensor]
                                                        else:
                                                                dataFlowSpeed[str(intervalIndex)].append(speedSensor)
                except IndexError:
                        pass

        collectedOccupancy = zeros(timeInterval-1)        
        for timeID in dataTime:
                if int(timeID) >= deltaT:
                        timeIndex = int(int(timeID)/deltaT)-1
                        try:
                                for j in range(len(dataTime[timeID][:,0])):
                                        if (sensorLocation - sensorLength/2 - 7) < dataTime[timeID][j,1] < (sensorLocation + sensorLength/2 + 7):
                                                collectedOccupancy[timeIndex] = collectedOccupancy[timeIndex] + 1

                        except IndexError:                                
                                if (sensorLocation - sensorLength/2 - 7) < dataTime[timeID][1] < (sensorLocation + sensorLength/2 + 7):
                                        collectedOccupancy[timeIndex] = collectedOccupancy[timeIndex] + 1
        collectedOccupancy = collectedOccupancy/3/deltaT ## per lane
        for k in dataFlowSpeed:
                flow = len(dataFlowSpeed[k])*(simulationTime/deltaT)/3
                speed = average(dataFlowSpeed[k])
                collectedData[k,0]=collectedOccupancy[int(k)]*5280/(14+sensorLength)
                collectedData[k,1]=speed
                collectedData[k,2]=flow
        return collectedData

if __name__=='__main__':

        ## the directory to load the extracted keyData
        loadDirectory = '/Users/Ren/Dropbox/SourceCode/CORSIM_FD/ExtractedKeyData/'
        ## the directory to store the processed keyData
        saveDirectory = '/Users/Ren/Dropbox/SourceCode/CORSIM_FD/ProcessedKeyData/'       
                
                
        fileList = ['L3F7000S65','L3F5000S65','L3F3000S65','L3F7000S45','L3F7000S35','L3F7000S25','L3F7000S15','L3F7000S10','L3F7000S5','L3F7000S1']
        ##fileList=['L3F3000S65']

        fileNumber = len(fileList)
        fileCount = 0
        for fileName in fileList:
                fileCount = fileCount+1
                print 'This is file', fileCount, 'out of', fileNumber
                
                keyData=load(loadDirectory + fileName+'key.npy')
                deltaT = 2*60 # second, this is the time interval for data collection
                simulationTime = 3600
                sensorLocation = 1000 # feet
                sensorLength = 6 # feet

                dataTime = build_data_time(keyData)
                vehTrajectory = build_veh_trajectory(keyData)
                collectedData = collect_flow_speed_occupancy(vehTrajectory,dataTime, deltaT, simulationTime, sensorLocation, sensorLength)


                print 'file'+fileName
                save(saveDirectory + fileName+'ProcessedkeyData',collectedData)


