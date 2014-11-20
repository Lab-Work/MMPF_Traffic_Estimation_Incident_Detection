## This function extracts the useful data from CORSIM trajectory text file
## The extracted trajectory data is saved as keyData, in the following format: [TimeIndex, VehID, VehPosition, Speed]
## The keyData file will be later used for fundamental diagram calibration purpose


from numpy import *
from copy import deepcopy

if __name__=='__main__':

        ## the directory where the output text files from CORSIM are located
        loadDirectory = '/Users/Ren/Dropbox/SourceCode/CORSIM_FD/CORSIMOutPutTextFile/'
        ## the directory where the output from this function will be saved at
        saveDirectory = '/Users/Ren/Dropbox//SourceCode/CORSIM_FD/ExtractedKeyData/'


        ## the list of files that will be extracted
        fileList = ['L3F7000S65','L3F5000S65','L3F3000S65','L3F7000S45','L3F7000S35','L3F7000S25','L3F7000S15','L3F7000S10','L3F7000S5','L3F7000S1',]
        fileList = ['L3F3000S65']
        fileNumber = len(fileList)

        fileCount=0
        for fileName in fileList:
                fileCount = fileCount+1
                print 'This is file', fileCount, 'out of', fileNumber

                keyData = zeros((10000000,4))
                lineNumber = 0
                index = 0
                f = open(loadDirectory + fileName+'.txt','r')
                for line in f:
                        lineNumber = lineNumber + 1
                        if mod(lineNumber,1000) == 0:
                                print 'This is iteration', lineNumber
                        nodeCount = 0
                        initialList = []
                        myFocus = line.split()
                        myFocus = [int(i) for i in myFocus]
                        if len(myFocus) < 10:
                                pass
                        else:
                                for i in range(len(myFocus)):
                                        if i == 0:
                                                initialList.append(myFocus[i+3])
                                        if (i-17) >= 0 and mod((i-17),19) == 0 and myFocus[i+1] == 3001 and (myFocus[i+1] != myFocus[i+2]):
                                                j=i
                                                break
                                for k in range(j,len(myFocus)):
                                        if (k-j-17) >= 0 and mod((k-j-17),19) == 0 and (myFocus[k+1] != 3003):
                                                currentList = list(initialList)
                                                currentList.append(myFocus[k+1])
                                                currentList.append(myFocus[k+8])
                                                currentList.append(myFocus[k+13]*3600/5280)
                                                keyData[index] = array(currentList)
                                                index = index+1
                                        elif (k-j-17) >= 0 and mod((k-j-17),19) == 0 and (myFocus[k+1] == 3003):
                                                break
                for i in reversed(range(len(keyData[:,0]))):
                        if keyData[i,3] != 0:
                                keyData = keyData[0:i+1,:].copy()
                                break
                save(saveDirectory + fileName + 'Key',keyData)
