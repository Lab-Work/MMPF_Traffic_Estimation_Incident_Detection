## This function extracts the useful data from CORSIM trajectory text file
## The extracted trajectory data is saved as keyData, in the following format: [TimeIndex, VehID, VehPosition, Speed, Lane]
## The keyData file will be later used for measurement generation

from numpy import *
from copy import deepcopy
import sys
import os

if __name__== '__main__':

        ## the directory where the output text files from CORSIM are located
        loadDirectory = os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/CORSIMOutPutTextFile/'
        ## the directory where the output from this function will be saved at
        saveDirectory = os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/ExtractedKeyData/'


        ## the list of files that will be extracted
#        fileList = ['L3B1F5000S65R50', 'L3B1F5000S65R50', 'L3B1F4000S65R50', 'L3B1F3000S65R50', 'L3B1F2000S65R50', 'L3B1F1000S65R50']
        fileList = ['L3B1F6000S65R50']
        fileNumber = len(fileList)

        fileCount = 0
        for fileName in fileList:
                fileCount = fileCount+1
                print 'This is file', fileCount, 'out of', fileNumber


                keyData = zeros((10000000,5))
                lineNumber = 0
                index = 0
                f=open(loadDirectory+fileName+'.txt','r')
                for line in f:
                        lineNumber = lineNumber + 1
                        if mod(lineNumber,1000) == 0:
                                print 'This is iteration', lineNumber
                        initialList=[]
                        myFocus = line.split()
                        myFocus = [int(i) for i in myFocus]
                        if len(myFocus)<10:
                                pass
                        else:
                                for i in range(int(len(myFocus))):
                                        if i == 3:
                                                initialList.append(myFocus[i])
                                        elif mod((i-18),19) == 0:
                                                if i+10>len(myFocus):
                                                        pass
                                                else:
                                                        currentList = list(initialList)
                                                        currentList.append(myFocus[i])
                                                        currentList.append(myFocus[i+7])
                                                        currentList.append(myFocus[i+12]*3600/5280.0)
                                                        currentList.append(myFocus[i+6])
                                                        keyData[index,:] = array(currentList)
                                                        index = index+1

                                
                for i in reversed(range(len(keyData[:,0]))):
                        if keyData[i,3] != 0:
                                keyData = keyData[0:i+1,:].copy()
                                break


                save(saveDirectory + fileName+'KeyTra',keyData)
