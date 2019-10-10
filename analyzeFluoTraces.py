import numpy as np
import matplotlib.pyplot as plt
import cv2
import pickle
import pdb
from ScanImageTiffReader import ScanImageTiffReader

############################################################
def readTimeStampOfRecording(tiffFile, nFrame):
    desc = ScanImageTiffReader(tiffFile).description(nFrame)
    keyWordIdx = desc.find('epoch')
    dateString = re.split('\[|\]', desc[keyWordIdx:])
    dateIdv = dateString[1].split()
    # print(dateIdv)
    unixStartTime = int(datetime.datetime(int(dateIdv[0]), int(dateIdv[1]), int(dateIdv[2]), int(dateIdv[3]), int(dateIdv[4]), int(float(dateIdv[5]))).strftime('%s'))
    #
    keyWordIdx = desc.find('frameTimestamps_sec')
    splitString = re.split('=|\n', desc[keyWordIdx:])
    frameTimestamps = float(splitString[1])

    keyWordIdx = desc.find('acqTriggerTimestamps_sec')
    splitString = re.split('=|\n', desc[keyWordIdx:])
    acqTriggerTimestamps = float(splitString[1])

    keyWordIdx = desc.find('frameNumberAcquisition')
    splitString = re.split('=|\n', desc[keyWordIdx:])
    frameNumberAcquisition = int(splitString[1])

    keyWordIdx = desc.find('acquisitionNumbers')
    splitString = re.split('=|\n', desc[keyWordIdx:])
    acquisitionNumbers = int(splitString[1])

    unixFrameTime = unixStartTime + frameTimestamps
    # print(tiffFile,unixTime)
    return ([frameNumberAcquisition, acquisitionNumbers, unixStartTime, unixFrameTime, frameTimestamps, acqTriggerTimestamps])

def extractAndSaveCaTimeStamps(saveDir,tiffPaths):
    timeStamps = []
    for i in range(len(tiffPaths)):
        data = ScanImageTiffReader(tiffPaths[i]).data()
        fN = np.shape(data)[0]
        #frameNumbers.append(fN)
        for n in range(fN):
            timeStamps.append(readTimeStampOfRecording(tiffPaths[i],n))

    timeStampsA = np.asarray(timeStamps)
    np.save(saveDir+'/suite2p/plane0/timeStamps.npy',timeStampsA)

#########################################################
# set parameters

dataOutDir = 'dataOutput/'
figOutDir = 'figureOutput/'

animalID = 'animal#1'
baseDir = '/media/labo_rw/JINmGluR/'
beforeDrugDir = '2019.08.15_002_animal#1/animal#1_00000-00013/'
afterDrugDir = '2019.08.15_002_animal#1/animal#1_00015-20_23-28/'
suite2pDir = 'suite2p/plane0/'

startBD = 0
endBD =13
nImgsBD = range(startBD,endBD+1)
nImgsAD = [15,16,17,18,19,20,23,24,25,26,27,28]


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

FBD = np.load(baseDir + beforeDrugDir + suite2pDir + 'F.npy')
FneuBD = np.load(baseDir + beforeDrugDir + suite2pDir + 'Fneu.npy')
FAD = np.load(baseDir + afterDrugDir + suite2pDir + 'F.npy')
FneuAD = np.load(baseDir + afterDrugDir + suite2pDir + 'Fneu.npy')

###########################################################
# read ROI matching results
intersectionROIs = pickle.load(open( dataOutDir + 'ROIintersections_%s.p' % animalID, 'rb'))

# read time stamps
timeStampsBD = []
for n in nImgsBD:
    tS = np.load(baseDir + beforeDrugDir + 'rawImages/%s_%05d_timeStamps.npy' % (animalID,n))
    timeStampsBD.append([n,tS])

timeStampsAD = []
for n in nImgsAD:
    tS = np.load(baseDir + afterDrugDir + 'rawImages/%s_%05d_timeStamps.npy' % (animalID,n))
    timeStampsAD.append([n,tS])

#pdb.set_trace()

##################################################################
# Show final results

fig_width = 12 # width in inches
fig_height = 8  # height in inches
fig_size =  [fig_width,fig_height]
params = {'axes.labelsize': 14,
          'axes.titlesize': 13,
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
# create figure instance
fig = plt.figure()

# define sub-panel grid and possibly width and height ratios
gs = gridspec.GridSpec(2, 2,
                       width_ratios=[1,1.2],
                       height_ratios=[1,1]
                       )
# define vertical and horizontal spacing between panels
gs.update(wspace=0.3,hspace=0.4)

# possibly change outer margins of the figure
#plt.subplots_adjust(left=0.14, right=0.92, top=0.92, bottom=0.18)

gssub0 = gridspec.GridSpecFromSubplotSpec(1, len(nImgsBD) , subplot_spec=gs[0],hspace=0.4)

for n in range(nImgsBD):
    ax0 = plt.subplot(gssub0[i])
    for i in range(len(intersectionROIs)):
        ax0.plot()


for i in range(len(intersectionROIs)):
    ax0.plot()


ax0 = plt.subplot(gs[1]) #####################
ax0.set_title('after drug')

for i in range(len(intersectionROIs)):



#plt.savefig(figOutDir+'FluorescenceTraces_%s.pdf' % animalID)
plt.show()