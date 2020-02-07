import numpy as np
import matplotlib.pyplot as plt
import cv2
import pickle
import pdb
from animalSettings import animalSettings
import sys


try:
    anim = sys.argv[1]
except IndexError:
    anim = 'animal#1_2'
else:
    pass

aS = animalSettings(anim)
aS.loadSettings()

print('animal : ',aS.animalID)

#########################################################
# set parameters

dataOutDir = 'dataOutput/'
figOutDir = 'figureOutput/'

redoIntersections = True


##################################################################
# extract and save time-stamps
# read time stamps
timeStampsBD = []
wheelBD = []
for i in range(len(aS.nImgsBD)):
    tS = np.load(aS.baseDir + aS.beforeDrugDir + 'rawImages/%s_%05d_timeStamps.npy' % (aS.animalID,aS.nImgsBD[i][0]))
    timeStampsBD.append([aS.nImgsBD[i][0],tS])
    wa = np.load(aS.baseDir + aS.beforeDrugDir + 'rawImages/walkingActivity_%s_rec%s.p' % (aS.beforeDrugDir[:14],aS.nImgsBD[i][1]))
    wheelBD.append([aS.nImgsBD[i][0],wa])

timeStampsAD = []
wheelAD = []
for i in range(len(aS.nImgsAD)):
    tS = np.load(aS.baseDir + aS.afterDrugDir + 'rawImages/%s_%05d_timeStamps.npy' % (aS.animalID,aS.nImgsAD[i][0]))
    #print(i,len(tS),aS.baseDir + aS.afterDrugDir + 'rawImages/%s_%05d_timeStamps.npy' % (aS.animalID,aS.nImgsAD[i][0]))
    timeStampsAD.append([aS.nImgsAD[i][0],tS])
    wa = np.load(aS.baseDir + aS.afterDrugDir + 'rawImages/walkingActivity_%s_rec%s.p' % (aS.afterDrugDir[:14], aS.nImgsAD[i][1]))
    wheelAD.append([aS.nImgsAD[i][0], wa])

timeStamps820 = []
wheel820 = []
for i in range(len(aS.nImgs820)):
    tS = np.load(aS.baseDir + aS.eightTwentyDir + 'rawImages/%s_%05d_timeStamps.npy' % (aS.animalID,aS.nImgs820[i][0]))
    #print(i,len(tS),aS.baseDir + aS.afterDrugDir + 'rawImages/%s_%05d_timeStamps.npy' % (aS.animalID,aS.nImgsAD[i][0]))
    timeStamps820.append([aS.nImgs820[i][0],tS])
    wa = np.load(aS.baseDir + aS.eightTwentyDir + 'rawImages/walkingActivity_%s_rec%s.p' % (aS.eightTwentyDir[:14], aS.nImgs820[i][1]))
    wheel820.append([aS.nImgs820[i][0], wa])


pickle.dump(timeStampsBD, open( dataOutDir + 'timeStampsBeforeDrug_%s.p' % aS.animalID, 'wb' ) )
pickle.dump(timeStampsAD, open( dataOutDir + 'timeStampsAfterDrug_%s.p' % aS.animalID, 'wb' ) )
pickle.dump(timeStamps820, open( dataOutDir + 'timeStamps820_%s.p' % aS.animalID, 'wb' ) )

pickle.dump(wheelBD, open( dataOutDir + 'wheelActivityBeforeDrug_%s.p' % aS.animalID, 'wb' ) )
pickle.dump(wheelAD, open( dataOutDir + 'wheelActivityAfterDrug_%s.p' % aS.animalID, 'wb' ) )
pickle.dump(wheel820, open( dataOutDir + 'wheelActivity820_%s.p' % aS.animalID, 'wb' ) )