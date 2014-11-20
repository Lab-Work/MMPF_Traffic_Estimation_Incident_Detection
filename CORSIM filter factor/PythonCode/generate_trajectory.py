## This file generates the true density and speed evolution of CORSIM simulations

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from numpy import *
from random import choice
from random import sample
from copy import deepcopy
import matplotlib.animation as animation


def build_data_time(keyData):
## construct the data time dictionary
## key: TimeIndex
## value: [VehID, VehPosition, Speed,laneID]
        dataTime = dict()
        for i in range(len(keyData[:,0])):
                if str(int(keyData[i,0])) not in dataTime:
                        dataTime[str(int(keyData[i,0]))] = keyData[i,[1,2,3,4]]
                else:
                        dataTime[str(int(keyData[i,0]))] = vstack((dataTime[str(int(keyData[i,0]))], keyData[i,[1,2,3,4]]))
        return dataTime

if __name__=='__main__':

        ## the directory to load the extracted keyData
        loadDirectory = '/Users/Ren/Dropbox/SourceCode/CORSIM filter factor/ExtractedKeyData/'
        ## the directory to store the true state
        saveDirectory = '/Users/Ren/Dropbox/SourceCode/CORSIM filter factor/TrueState/'       

        fileName = 'L3B1F6000S65R50'
        keyData = load(loadDirectory+fileName+'keyTra.npy')
        dataTime = build_data_time(keyData)

        plotlays = [0,0,0]
        plotcols = ['blue','blue','blue']
        label = ['free flow','medium flow','congested']

        fig = plt.figure()
        ax = plt.axes(xlim=(0, 4.0*5280), ylim=(0.5, 3.5))
#        ax = plt.axes(xlim=(8548, 9548), ylim=(0, 4))

        timetext = ax.text(1000,3.2,'')
        incidenttext = ax.text(7800,1.2,'')
        sensortext1 = ax.text(0.4*5280,0.3,'Loop detector')
        sensortext2 = ax.text(3.6*5280,0.3,'Loop detector')

        lines = []


        for index,lay in enumerate(plotlays):
                if index == 0:
                        lobj = ax.plot([],[],'ro',lw=2,color=plotcols[index],label=label[index])[0]
                        lines.append(lobj)
                elif index == 1:
                        lobj = ax.plot([],[],'ro',lw=2,color=plotcols[index],label=label[index])[0]
                        lines.append(lobj)
                elif index ==2:
                        lobj = ax.plot([],[],'ro',lw=2,color=plotcols[index],label=label[index])[0]
                        lines.append(lobj)
        ax.legend()
        plt.rc('xtick',labelsize=20)
        plt.rc('ytick',labelsize=20)
        plt.xlabel('space',fontsize=20)
        plt.ylabel('lane ID',fontsize=20)
        plt.yticks([1.0, 2.0,3.0])
        plt.xticks([])
        plt.title('vehicle trajectory',fontsize=20)


        def init():
            for line in lines:
                line.set_data([],[])
            return lines

        def animate(i):
            timetext.set_text('Seconds: '+str(i))
            if i>=1200 and i<=2400:
                    incidenttext.set_text('Incident occurs')
            elif i>2400:
                    incidenttext.set_text('Incident cleared')
        

            try:
                    x = dataTime[str(i)][:,1]
                    y = dataTime[str(i)][:,3]
            except IndexError:
                    x = dataTime[str(i)][1]
                    y = dataTime[str(i)][3]

            for lnum,line in enumerate(lines):
                line.set_data(x,y)

            return lines, timetext

        anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=3600, interval=100, blit=False, repeat = False)

        plt.show()
        anim.save('trajectoryShort'+'.mp4',fps=100)


