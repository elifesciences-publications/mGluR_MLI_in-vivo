import numpy as np
import matplotlib.pyplot as plt
import cv2
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

##########################################################
# read data and determine principal parameters
statBD = np.load(aS.baseDir + aS.beforeDrugDir + aS.suite2pDir + 'stat.npy')
opsBD = np.load(aS.baseDir + aS.beforeDrugDir + aS.suite2pDir + 'ops.npy').item()

statAD = np.load(aS.baseDir + aS.afterDrugDir + aS.suite2pDir + 'stat.npy')
opsAD = np.load(aS.baseDir + aS.afterDrugDir + aS.suite2pDir + 'ops.npy').item()


emptyImBD = np.zeros((opsBD['Ly'], opsBD['Lx']))
emptyImAD = np.zeros((opsAD['Ly'], opsAD['Lx']))

# Read the enhanced ROI images to be aligned
imBD =  opsBD['meanImgE'][aS.cutOffX:-aS.cutOffX,aS.cutOffY:-(aS.cutOffY+1)]
imAD =  opsAD['meanImgE'][aS.cutOffX:-aS.cutOffX,aS.cutOffY:-(aS.cutOffY+1)]

###########################################################
# perform image alignment

# Convert images to grayscale
#im1_gray = cv2.cvtColor(im1, cv2.COLOR_BGR2GRAY)
#im2_gray = cv2.cvtColor(im2, cv2.COLOR_BGR2GRAY)

# Find size of image1
sz = imBD.shape

# Define the motion model
warp_mode = cv2.MOTION_TRANSLATION # MOTION_EUCLIDEAN

# Define 2x3 or 3x3 matrices and initialize the matrix to identity
if warp_mode == cv2.MOTION_HOMOGRAPHY:
    warp_matrix = np.eye(3, 3, dtype=np.float32)
else:
    warp_matrix = np.eye(2, 3, dtype=np.float32)

warp_matrix[0,2] = aS.xOffset
warp_matrix[1,2] = aS.yOffset


# Specify the number of iterations.
number_of_iterations = 5000;

# Specify the threshold of the increment
# in the correlation coefficient between two iterations
termination_eps = 1e-10;

# Define termination criteria
criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, number_of_iterations, termination_eps)

# Run the ECC algorithm. The results are stored in warp_matrix.
(cc, warp_matrix) = cv2.findTransformECC(imBD, imAD, warp_matrix, warp_mode, criteria)

if warp_mode == cv2.MOTION_HOMOGRAPHY:
    # Use warpPerspective for Homography
    imAD_aligned = cv2.warpPerspective(imAD, warp_matrix, (sz[1], sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP)
else:
    # Use warpAffine for Translation, Euclidean and Affine
    imAD_aligned = cv2.warpAffine(imAD, warp_matrix, (sz[1], sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP);

print('result of image alignment-> warp-matrix and correlation coefficient : ', warp_matrix, cc)

np.save(dataOutDir+'warp_matrix_%s.npy' % aS.animalID, warp_matrix)

# for n in range(0,ncells):
#     ypix = stat[n]['ypix'][~stat[n]['overlap']]
#     xpix = stat[n]['xpix'][~stat[n]['overlap']]
#     im[ypix,xpix] = n+1

##################################################################
# Show final results
fig = plt.figure(figsize=(10,10))

plt.figtext(0.1, 0.95, '%s ' % (aS.animalID),clip_on=False, color='black', size=14)

ax0 = fig.add_subplot(2,2,1)
ax0.set_title('before drug')
ax0.imshow(imBD)

ax0 = fig.add_subplot(2,2,2)
ax0.set_title('after drug')
ax0.imshow(imAD)

ax0 = fig.add_subplot(2,2,3)
ax0.set_title('overlay of original images')
overlayBefore = cv2.addWeighted(imBD, 1, imAD, 1, 0)
ax0.imshow(overlayBefore)

ax0 = fig.add_subplot(2,2,4)
ax0.set_title('overlay after alignement of images')
overlayAfter = cv2.addWeighted(imBD, 1, imAD_aligned, 1, 0)
ax0.imshow(overlayAfter)


plt.savefig(figOutDir+'ImageAlignment_%s.pdf' % aS.animalID)
#plt.savefig(figOutDir+'ImageAlignment_%s.png' % aS.animalID)
#plt.show()