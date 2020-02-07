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

stat820 = np.load(aS.baseDir + aS.eightTwentyDir + aS.suite2pDir + 'stat.npy')
ops820 = np.load(aS.baseDir + aS.eightTwentyDir + aS.suite2pDir + 'ops.npy').item()
iscell820 = np.load(aS.baseDir + aS.eightTwentyDir + aS.suite2pDir +'iscell.npy')
ncells820 = len(iscell820)

emptyImBD = np.zeros((opsBD['Ly'], opsBD['Lx']))
emptyImAD = np.zeros((opsAD['Ly'], opsAD['Lx']))
emptyIm820 = np.zeros((ops820['Ly'], ops820['Lx']))

# Read the enhanced ROI images to be aligned
imBD =  opsBD['meanImgE'][aS.cutOffX:-aS.cutOffX,aS.cutOffY:-(aS.cutOffY+1)]
imAD =  opsAD['meanImgE'][aS.cutOffX:-aS.cutOffX,aS.cutOffY:-(aS.cutOffY+1)]
im820 =  ops820['meanImgE'][aS.cutOffX:-aS.cutOffX,aS.cutOffY:-(aS.cutOffY+1)]

warp_matrix1 = np.load(dataOutDir+'warp_matrix_%s.npy' % aS.animalID)
warp_matrix2 = np.load(dataOutDir+'warp_matrix820_%s.npy' % aS.animalID)

###########################################################
# read and match ROIs

imMaskBD =  np.zeros((opsBD['Ly'], opsBD['Lx']))
imMaskAD = np.zeros((opsAD['Ly'], opsAD['Lx']))
imMask820 = np.zeros((ops820['Ly'], ops820['Lx']))

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
                    xpixADPrime = xpixAD - np.rint(warp_matrix1[0,2])
                    ypixADPrime = ypixAD - np.rint(warp_matrix1[1,2])
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

    #############
    intersectionROIs_820 = []
    for n in range(0, ncellsBD):
        imMaskBD[:] = 0
        if iscellBD[n][0] == 1:
            ypixBD = statBD[n]['ypix']
            xpixBD = statBD[n]['xpix']
            imMaskBD[ypixBD, xpixBD] = 1
            for m in range(0, ncells820):
                imMask820[:] = 0
                if iscell820[m][0] == 1:
                    ypix820 = stat820[m]['ypix']
                    xpix820 = stat820[m]['xpix']
                    # perform eucledian transform
                    xpix820Prime = xpix820 - np.rint(warp_matrix2[0, 2])
                    ypix820Prime = ypix820 - np.rint(warp_matrix2[1, 2])
                    xpix820Prime = np.array(xpix820Prime, dtype=int)
                    ypix820Prime = np.array(ypix820Prime, dtype=int)
                    # make sure pixels remain within
                    xpix820Prime2 = xpix820Prime[(xpix820Prime < ops820['Lx']) & (ypix820Prime < ops820['Ly'])]
                    ypix820Prime2 = ypix820Prime[(xpix820Prime < ops820['Lx']) & (ypix820Prime < ops820['Ly'])]
                    imMask820[ypix820Prime2, xpix820Prime2] = 1
                    intersection = np.sum(np.logical_and(imMaskBD, imMask820))
                    eitherOr = np.sum(np.logical_or(imMaskBD, imMask820))
                    if intersection > 0:
                        print(n, m, intersection, eitherOr, intersection / eitherOr)
                        intersectionROIs_820.append([n, m, xpixBD, ypixBD, xpix820Prime2, ypix820Prime2, intersection, eitherOr, intersection / eitherOr])

    pickle.dump(intersectionROIs_820, open(dataOutDir + 'ROIintersections820_%s.p' % aS.animalID, 'wb'))

else:
    intersectionROIs = pickle.load(open( dataOutDir + 'ROIintersections_%s.p' % aS.animalID, 'rb'))
    intersectionROIs_820 = pickle.load(open( dataOutDir + 'ROIintersections820_%s.p' % aS.animalID, 'rb'))



##################################################################
# Show final results
fig = plt.figure(figsize=(15,15)) ########################

plt.figtext(0.1, 0.95, '%s ' % (aS.animalID),clip_on=False, color='black', size=14)

ax0 = fig.add_subplot(4,3,1) #############################
ax0.set_title('before drug')
ax0.imshow(imBD)

ax0 = fig.add_subplot(4,3,2) #############################
ax0.set_title('after drug')
ax0.imshow(imAD)

ax0 = fig.add_subplot(4,3,3) #############################
ax0.set_title('820 nm run')
ax0.imshow(im820)

ax0 = fig.add_subplot(4,3,4) #############################
ax0.set_title('ROIs before drug')
imBD = np.zeros((opsBD['Ly'], opsBD['Lx']))

for n in range(0,ncellsBD):
    if iscellBD[n][0]==1:
        ypixBD = statBD[n]['ypix']
        xpixBD = statBD[n]['xpix']
        imBD[ypixBD,xpixBD] = n+1

ax0.imshow(imBD,cmap='gist_ncar')

ax0 = fig.add_subplot(4,3,5) #############################
ax0.set_title('ROIs after drug')
imAD = np.zeros((opsAD['Ly'], opsAD['Lx']))

for n in range(0,ncellsAD):
    if iscellAD[n][0]==1:
        ypixAD = statAD[n]['ypix']
        xpixAD = statAD[n]['xpix']
        imAD[ypixAD,xpixAD] = n+1

ax0.imshow(imAD,cmap='gist_ncar')

ax0 = fig.add_subplot(4,3,6) #############################
ax0.set_title('ROIs at 820')
im820 = np.zeros((ops820['Ly'], ops820['Lx']))

for n in range(0,ncells820):
    if iscell820[n][0]==1:
        ypix820 = stat820[n]['ypix']
        xpix820 = stat820[n]['xpix']
        im820[ypix820,xpix820] = n+1

ax0.imshow(im820,cmap='gist_ncar')


ax0 = fig.add_subplot(4,3,8)  #############################
ax0.set_title('overlapping ROIs BD-AD')
imBD = np.zeros((opsBD['Ly'], opsBD['Lx']))
imAD = np.zeros((opsAD['Ly'], opsAD['Lx']))

for n in range(0,len(intersectionROIs)):

    ypixBD = intersectionROIs[n][3]
    xpixBD = intersectionROIs[n][2]
    ypixAD = intersectionROIs[n][5]
    xpixAD = intersectionROIs[n][4]
    imBD[ypixBD,xpixBD] = 1
    imAD[ypixAD,xpixAD] = 2

overlayBothROIs1 = cv2.addWeighted(imBD, 1, imAD, 1, 0)
ax0.imshow(overlayBothROIs1)

ax0 = fig.add_subplot(4,3,9)  #############################
ax0.set_title('overlapping ROIs BD-820')
imBD = np.zeros((opsBD['Ly'], opsBD['Lx']))
im820 = np.zeros((ops820['Ly'], ops820['Lx']))

for n in range(0,len(intersectionROIs_820)):

    ypixBD = intersectionROIs_820[n][3]
    xpixBD = intersectionROIs_820[n][2]
    ypix820 = intersectionROIs_820[n][5]
    xpix820 = intersectionROIs_820[n][4]
    imBD[ypixBD,xpixBD] = 1
    im820[ypix820,xpix820] = 2

overlayBothROIs2 = cv2.addWeighted(imBD, 1, im820, 1, 0)
ax0.imshow(overlayBothROIs2)

ax0 = fig.add_subplot(4,3,11)  #############################
ax0.set_title('fraction of ROI overlap BD-AD')
interFractions1 = []
for n in range(0,len(intersectionROIs)):
    interFractions1.append(intersectionROIs[n][8])

ax0.hist(interFractions1,bins=15)

ax0 = fig.add_subplot(4,3,12)  #############################
ax0.set_title('fraction of ROI overlap BD-820')
interFractions2 = []
for n in range(0,len(intersectionROIs_820)):
    interFractions2.append(intersectionROIs_820[n][8])

ax0.hist(interFractions2,bins=15)



plt.savefig(figOutDir+'ROImatching_%s.pdf' % aS.animalID)
#plt.savefig(figOutDir+'ROImatching_%s.png' % aS.animalID)
#plt.show()