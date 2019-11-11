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
    anim = 'animal#2'
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

##########################################################
# read data and determine principal parameters
statBD = np.load(aS.baseDir + aS.beforeDrugDir + aS.suite2pDir + 'stat.npy')
opsBD = np.load(aS.baseDir + aS.beforeDrugDir + aS.suite2pDir + 'ops.npy').item()
iscellBD = np.load(aS.baseDir + aS.beforeDrugDir + aS.suite2pDir +'iscell.npy')
ncellsBD = len(iscellBD)

statAD = np.load(aS.baseDir + aS.afterDrugDir + aS.suite2pDir + 'stat.npy')
opsAD = np.load(aS.baseDir + aS.afterDrugDir + aS.suite2pDir + 'ops.npy').item()
iscellAD = np.load(aS.baseDir + aS.afterDrugDir + aS.suite2pDir +'iscell.npy')
ncellsAD = len(iscellAD)

emptyImBD = np.zeros((opsBD['Ly'], opsBD['Lx']))
emptyImAD = np.zeros((opsAD['Ly'], opsAD['Lx']))

# Read the enhanced ROI images to be aligned
imBD =  opsBD['meanImgE'][aS.cutOffX:-aS.cutOffX,aS.cutOffY:-(aS.cutOffY+1)]
imAD =  opsAD['meanImgE'][aS.cutOffX:-aS.cutOffX,aS.cutOffY:-(aS.cutOffY+1)]

warp_matrix = np.load(dataOutDir+'warp_matrix_%s.npy' % aS.animalID)

###########################################################
# read and match ROIs

imMaskBD =  np.zeros((opsBD['Ly'], opsBD['Lx']))
imMaskAD = np.zeros((opsAD['Ly'], opsAD['Lx']))

if redoIntersections:
    intersectionROIs = []
    for n in range(0,ncellsBD):
        imMaskBD[:] = 0
        if iscellBD[n][0]==1:
            ypixBD = statBD[n]['ypix']
            xpixBD = statBD[n]['xpix']
            imMaskBD[ypixBD,xpixBD] = 1
            for m in range(0,ncellsAD):
                imMaskAD[:] = 0
                if iscellAD[m][0]==1:
                    ypixAD = statAD[m]['ypix']
                    xpixAD = statAD[m]['xpix']
                    # perform eucledian transform
                    xpixADPrime = xpixAD - np.rint(warp_matrix[0,2])
                    ypixADPrime = ypixAD - np.rint(warp_matrix[1,2])
                    xpixADPrime = np.array(xpixADPrime,dtype=int)
                    ypixADPrime = np.array(ypixADPrime, dtype=int)
                    # make sure pixels remain within
                    xpixADPrime2 = xpixADPrime[(xpixADPrime<opsAD['Lx'])&(ypixADPrime<opsAD['Ly'])]
                    ypixADPrime2 = ypixADPrime[(xpixADPrime<opsAD['Lx'])&(ypixADPrime<opsAD['Ly'])]
                    imMaskAD[ypixADPrime2,xpixADPrime2] = 1
                    intersection = np.sum(np.logical_and(imMaskBD,imMaskAD))
                    eitherOr = np.sum(np.logical_or(imMaskBD,imMaskAD))
                    if intersection>0:
                        print(n,m,intersection,eitherOr,intersection/eitherOr)
                        intersectionROIs.append([n,m,xpixBD,ypixBD,xpixADPrime2,ypixADPrime2,intersection,eitherOr,intersection/eitherOr])

    pickle.dump(intersectionROIs, open( dataOutDir + 'ROIintersections_%s.p' % aS.animalID, 'wb' ) )
else:
    intersectionROIs = pickle.load(open( dataOutDir + 'ROIintersections_%s.p' % aS.animalID, 'rb'))

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
    wa = np.load(aS.baseDir + aS.afterDrugDir + 'rawImages/walkingActivity_%s_rec%s.p' % (aS.beforeDrugDir[:14], aS.nImgsAD[i][1]))
    wheelAD.append([aS.nImgsAD[i][0], wa])

pickle.dump(timeStampsBD, open( dataOutDir + 'timeStampsBeforeDrug_%s.p' % aS.animalID, 'wb' ) )
pickle.dump(timeStampsAD, open( dataOutDir + 'timeStampsAfterDrug_%s.p' % aS.animalID, 'wb' ) )

pickle.dump(wheelBD, open( dataOutDir + 'wheelActivityBeforeDrug_%s.p' % aS.animalID, 'wb' ) )
pickle.dump(wheelAD, open( dataOutDir + 'wheelActivityAfterDrug_%s.p' % aS.animalID, 'wb' ) )

##################################################################
# Show final results
fig = plt.figure(figsize=(10,15)) ########################

plt.figtext(0.1, 0.95, '%s ' % (aS.animalID),clip_on=False, color='black', size=14)

ax0 = fig.add_subplot(3,2,1) #############################
ax0.set_title('before drug')
ax0.imshow(imBD)

ax0 = fig.add_subplot(3,2,2) #############################
ax0.set_title('after drug')
ax0.imshow(imAD)

ax0 = fig.add_subplot(3,2,3) #############################
ax0.set_title('ROIs before drug')
imBD = np.zeros((opsBD['Ly'], opsBD['Lx']))

for n in range(0,ncellsBD):
    if iscellBD[n][0]==1:
        ypixBD = statBD[n]['ypix']
        xpixBD = statBD[n]['xpix']
        imBD[ypixBD,xpixBD] = n+1

ax0.imshow(imBD,cmap='gist_ncar')

ax0 = fig.add_subplot(3,2,4) #############################
ax0.set_title('ROIs after drug')
imAD = np.zeros((opsAD['Ly'], opsAD['Lx']))

for n in range(0,ncellsAD):
    if iscellAD[n][0]==1:
        ypixAD = statAD[n]['ypix']
        xpixAD = statAD[n]['xpix']
        imAD[ypixAD,xpixAD] = n+1

ax0.imshow(imAD,cmap='gist_ncar')


ax0 = fig.add_subplot(3,2,5)  #############################
ax0.set_title('overlapping ROIs')
imBD = np.zeros((opsBD['Ly'], opsBD['Lx']))
imAD = np.zeros((opsAD['Ly'], opsAD['Lx']))

for n in range(0,len(intersectionROIs)):

    ypixBD = intersectionROIs[n][3]
    xpixBD = intersectionROIs[n][2]
    ypixAD = intersectionROIs[n][5]
    xpixAD = intersectionROIs[n][4]
    imBD[ypixBD,xpixBD] = 1
    imAD[ypixAD,xpixAD] = 2

overlayBothROIs = cv2.addWeighted(imBD, 1, imAD, 1, 0)
ax0.imshow(overlayBothROIs)

ax0 = fig.add_subplot(3,2,6)  #############################
ax0.set_title('fraction of ROI overlap')
interFractions = []
for n in range(0,len(intersectionROIs)):
    interFractions.append(intersectionROIs[n][8])

ax0.hist(interFractions,bins=15)

plt.savefig(figOutDir+'ROImatching_%s.pdf' % aS.animalID)
#plt.savefig(figOutDir+'ROImatching_%s.png' % aS.animalID)
#plt.show()