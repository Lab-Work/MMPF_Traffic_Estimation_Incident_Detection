import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D
from numpy import *
from random import choice
from copy import deepcopy
import timeit
import sys
import os



flow = 6000
Density = load(os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir))+'/CORSIM filter factor/Measurements/meaDensityL3B1F'+str(flow)+'S65R50.npy')
VL = 14.0
DL = 6.0
Occupancy = Density*(VL+DL)/5280.0/3

T_OCCDF = 0.27
T_OCCRDF = 0.55
T_OCCTDF = 0.0003

TimeStep=179
OccCollected=zeros((30,2))

j=0
for i in range(TimeStep):
        if mod(i,6)==0:
                OccCollected[j,0] = Occupancy[i,1]
                OccCollected[j,1] = Occupancy[i,-2]
                j=j+1
                
for k in range(1,30):
        OCCDF = OccCollected[k,0]-OccCollected[k,1]
        OCCRDF = OCCDF/OccCollected[k,0]
        OCCTDF = (OccCollected[k-1,1]-OccCollected[k,1])/OccCollected[k-1,1]
##        print 'time step',k
##        print 'OCCDF', OCCDF
##        print 'OCCRDF', OCCRDF
##        print 'OCCTDF', OCCTDF
        if OCCDF>T_OCCDF and OCCRDF>T_OCCRDF and OCCTDF>T_OCCTDF:
                print 'Incident at time step', k*6
        


