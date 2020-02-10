import numpy as np
import cv2
import pdb
from animalSettings import animalSettings
import sys
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


try:
    anim = sys.argv[1]
except IndexError:
    anim = 'animal#1'
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

stat820 = np.load(aS.baseDir + aS.eightTwentyDir + aS.suite2pDir + 'stat.npy')
ops820 = np.load(aS.baseDir + aS.eightTwentyDir + aS.suite2pDir + 'ops.npy').item()


emptyImBD = np.zeros((opsBD['Ly'], opsBD['Lx']))
emptyImAD = np.zeros((opsAD['Ly'], opsAD['Lx']))
emptyIm820 = np.zeros((ops820['Ly'], ops820['Lx']))

# Read the enhanced ROI images to be aligned
imBD =  opsBD['meanImgE'][aS.cutOffX:-aS.cutOffX,aS.cutOffY:-(aS.cutOffY+1)]
imAD =  opsAD['meanImgE'][aS.cutOffX:-aS.cutOffX,aS.cutOffY:-(aS.cutOffY+1)]
im820 = ops820['meanImgE'][aS.cutOffX:-aS.cutOffX,aS.cutOffY:-(aS.cutOffY+1)]

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
    warp_matrix1 = np.eye(3, 3, dtype=np.float32)
    warp_matrix820 = np.eye(3, 3, dtype=np.float32)
else:
    warp_matrix1 = np.eye(2, 3, dtype=np.float32)
    warp_matrix820 = np.eye(2, 3, dtype=np.float32)

warp_matrix1[0,2] = aS.xOffset
warp_matrix1[1,2] = aS.yOffset

warp_matrix820[0,2] = aS.xOffset820
warp_matrix820[1,2] = aS.yOffset820

#warp_matrix2[0,2] = 0.
#warp_matrix2[1,2] = 0.

# Specify the number of iterations.
number_of_iterations = 5000;

# Specify the threshold of the increment
# in the correlation coefficient between two iterations
termination_eps = 1e-10;

# Define termination criteria
criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, number_of_iterations, termination_eps)

# Run the ECC algorithm. The results are stored in warp_matrix.
(cc1, warp_matrix1Ret) = cv2.findTransformECC(imBD, imAD, warp_matrix1, warp_mode, criteria)
(cc2,warp_matrix2Ret) = cv2.findTransformECC(imBD, im820, warp_matrix820, warp_mode, criteria)

if warp_mode == cv2.MOTION_HOMOGRAPHY:
    # Use warpPerspective for Homography
    imAD_aligned = cv2.warpPerspective(imAD, warp_matrix, (sz[1], sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP)
else:
    # Use warpAffine for Translation, Euclidean and Affine
    imAD_aligned = cv2.warpAffine(imAD, warp_matrix1Ret, (sz[1], sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP);
    im820_aligned = cv2.warpAffine(im820, warp_matrix2Ret, (sz[1], sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP);

print('result of image alignment-> warp-matrix BD-AD and correlation coefficient : ', warp_matrix1Ret, cc1)
print('result of image alignment-> warp-matrix BD-820 and correlation coefficient : ', warp_matrix2Ret, cc2)

np.save(dataOutDir+'warp_matrix_%s.npy' % aS.animalID, warp_matrix1Ret)
np.save(dataOutDir+'warp_matrix820_%s.npy' % aS.animalID, warp_matrix2Ret)

# for n in range(0,ncells):
#     ypix = stat[n]['ypix'][~stat[n]['overlap']]
#     xpix = stat[n]['xpix'][~stat[n]['overlap']]
#     im[ypix,xpix] = n+1

##################################################################
# Show final results
# figure #################################
fig_width = 10 # width in inches
fig_height = 10  # height in inches
fig_size =  [fig_width,fig_height]
params = {'axes.labelsize': 11,
          'axes.titlesize': 11,
          'font.size': 11,
          'xtick.labelsize': 11,
          'ytick.labelsize': 11,
          'figure.figsize': fig_size,
          'savefig.dpi' : 600,
          'axes.linewidth' : 1.3,
          'ytick.major.size' : 4,      # major tick size in points
          'xtick.major.size' : 4      # major tick size in points
          #'edgecolor' : None
          #'xtick.major.size' : 2,
          #'ytick.major.size' : 2,
          }
rcParams.update(params)

# set sans-serif font to Arial
rcParams['font.sans-serif'] = 'Arial'

# create figure instance
fig = plt.figure()


# define sub-panel grid and possibly width and height ratios
gs = gridspec.GridSpec(3, 3#,
                       #width_ratios=[1.2,1]
                       #height_ratios=[1,1]
                       )

# define vertical and horizontal spacing between panels
gs.update(wspace=0.3,hspace=0.3)

# possibly change outer margins of the figure
plt.subplots_adjust(left=0.05, right=0.95, top=0.92, bottom=0.06)

# sub-panel enumerations
#plt.figtext(0.06, 0.92, 'A',clip_on=False,color='black', weight='bold',size=22)


# first sub-plot #######################################################
# gssub = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0],hspace=0.2)
#ax0 = plt.subplot(gssub[0])

#fig = plt.figure(figsize=(10,10))

plt.figtext(0.1, 0.95, '%s ' % (aS.animalID),clip_on=False, color='black', size=14)

ax0 = plt.subplot(gs[0])
ax0.set_title('before drug')
ax0.imshow(imBD)

ax0 = plt.subplot(gs[1])
ax0.set_title('after drug')
ax0.imshow(imAD)

ax0 = plt.subplot(gs[2])
ax0.set_title('820 nm')
ax0.imshow(im820)

ax0 = plt.subplot(gs[4])
ax0.set_title('overlay of BD-AD images')
overlayBefore = cv2.addWeighted(imBD, 1, imAD, 1, 0)
ax0.imshow(overlayBefore)

ax0 = plt.subplot(gs[5])
ax0.set_title('overlay of BD-820 images')
overlayBefore = cv2.addWeighted(imBD, 1, im820, 1, 0)
ax0.imshow(overlayBefore)


ax0 = plt.subplot(gs[7])
ax0.set_title('overlay after alignement \nof BD-AD images',fontsize=10)
overlayAfter = cv2.addWeighted(imBD, 1, imAD_aligned, 1, 0)
ax0.imshow(overlayAfter)

ax0 = plt.subplot(gs[8])
ax0.set_title('overlay after alignement \nof BD-820 images',fontsize=10)
overlayAfter = cv2.addWeighted(imBD, 1, im820_aligned, 1, 0)
ax0.imshow(overlayAfter)


plt.savefig(figOutDir+'ImageAlignment_%s.pdf' % aS.animalID)
#plt.savefig(figOutDir+'ImageAlignment_%s.png' % aS.animalID)
#plt.show()