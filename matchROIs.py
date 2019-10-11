import numpy as np
import matplotlib.pyplot as plt
import cv2
import pickle
import pdb

#########################################################
# set parameters

dataOutDir = 'dataOutput/'
figOutDir = 'figureOutput/'

animalID = 'animal#1'
baseDir = '/media/labo_rw/JINmGluR/'
beforeDrugDir = '2019.08.15_002_animal#1/animal#1_00000-00013/'
afterDrugDir = '2019.08.15_002_animal#1/animal#1_00015-20_23-28/'
suite2pDir = 'suite2p/plane0/'

cutOffX = 8
cutOffY = 19

redoIntersections = False

##########################################################
# read data and determine principal parameters
statBD = np.load(baseDir + beforeDrugDir + suite2pDir + 'stat.npy')
opsBD = np.load(baseDir + beforeDrugDir + suite2pDir + 'ops.npy').item()
iscellBD = np.load(baseDir + beforeDrugDir + suite2pDir +'iscell.npy')
ncellsBD = len(iscellBD)

statAD = np.load(baseDir + afterDrugDir + suite2pDir + 'stat.npy')
opsAD = np.load(baseDir + afterDrugDir + suite2pDir + 'ops.npy').item()
iscellAD = np.load(baseDir + afterDrugDir + suite2pDir +'iscell.npy')
ncellsAD = len(iscellAD)

emptyImBD = np.zeros((opsBD['Ly'], opsBD['Lx']))
emptyImAD = np.zeros((opsAD['Ly'], opsAD['Lx']))

# Read the enhanced ROI images to be aligned
imBD =  opsBD['meanImgE'][cutOffX:-cutOffX,cutOffY:-(cutOffY+1)]
imAD =  opsAD['meanImgE'][cutOffX:-cutOffX,cutOffY:-(cutOffY+1)]

warp_matrix = np.load(dataOutDir+'warp_matrix_%s.npy' % animalID)

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

    pickle.dump(intersectionROIs, open( dataOutDir + 'ROIintersections_%s.p' % animalID, 'wb' ) )
else:
    intersectionROIs = pickle.load(open( dataOutDir + 'ROIintersections_%s.p' % animalID, 'rb'))

##################################################################
# extract and save time-stamps
# read time stamps
timeStampsBD = []
for n in nImgsBD:
    tS = np.load(baseDir + beforeDrugDir + 'rawImages/%s_%05d_timeStamps.npy' % (animalID,n))
    timeStampsBD.append([n,tS])

timeStampsAD = []
for n in nImgsAD:
    tS = np.load(baseDir + afterDrugDir + 'rawImages/%s_%05d_timeStamps.npy' % (animalID,n))
    timeStampsAD.append([n,tS])

pickle.dump(timeStampsBD, open( dataOutDir + 'timeStampsBeforeDrug_%s.p' % animalID, 'wb' ) )
pickle.dump(timeStampsAD, open( dataOutDir + 'timeStampsAfterDrug_%s.p' % animalID, 'wb' ) )

##################################################################
# Show final results
fig = plt.figure(figsize=(10,15)) ########################

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

plt.savefig(figOutDir+'ROImatching_%s.pdf' % animalID)
plt.show()