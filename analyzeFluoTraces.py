import numpy as np
import matplotlib.pyplot as plt
import cv2
import pickle
import pdb
from ScanImageTiffReader import ScanImageTiffReader
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
from scipy import integrate

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
#baseDir = '/media/labo_rw/JINmGluR/'
baseDir = '/home/mgraupe/'
beforeDrugDir = '2019.08.15_002_animal#1/animal#1_00000-00013/'
afterDrugDir = '2019.08.15_002_animal#1/animal#1_00015-20_23-28/'
suite2pDir = 'suite2p/plane0/'

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

#################################################################
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

#################################################################
# calculate dF/F0 and integral

minOverlap = 0.5
startN  = 0
baselineTime = 6.

analysisBD = []
integrals = []
for i in range(len(intersectionROIs)):
    if intersectionROIs[i][8]>minOverlap:
        startNBD = 0
        iE = []
        for n in range(len(nImgsBD)):
            tBD = timeStampsBD[n][1][:,1]
            FBDsingle = FBD[intersectionROIs[i][0]][startNBD:(startNBD+len(timeStampsBD[n][1][:,1]))] - 0.7*FneuBD[intersectionROIs[i][0]][startNBD:(startNBD+len(timeStampsBD[n][1][:,1]))]
            baselineMask = tBD < baselineTime
            F0BD = np.mean(FBDsingle[baselineMask])
            deltaBD = (FBDsingle-F0BD)/F0BD
            integralBD = integrate.simps(deltaBD, tBD)
            #analysisBD.append([n,i,intersectionROIs[i],tBD,FBD,F0BD,deltaBD,integralBD])
            iE.append(integralBD)
            startNBD += len(timeStampsBD[n][1])
        startNAD = 0
        for n in range(len(nImgsAD)):
            tAD = timeStampsAD[n][1][:,1]
            FADsingle = FAD[intersectionROIs[i][1]][startNAD:(startNAD+len(timeStampsAD[n][1][:,1]))] - 0.7*FneuAD[intersectionROIs[i][1]][startNAD:(startNAD+len(timeStampsAD[n][1][:,1]))]
            baselineMask = tAD < baselineTime
            F0AD = np.mean(FADsingle[baselineMask])
            deltaAD = (FADsingle-F0AD)/F0AD
            integralAD = integrate.simps(deltaAD, tAD)
            #analysisBD.append([n,i,intersectionROIs[i],tBD,FBD,F0BD,deltaBD,integralBD])
            iE.append(integralAD)
            startNAD += len(timeStampsAD[n][1])
        integrals.append(iE)

integrals = np.asarray(integrals)
#pdb.set_trace()
for i in range(len(integrals)):
    plt.plot([np.mean(integrals[i][:14]),np.mean(integrals[i][14:])])

plt.show()
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
gs = gridspec.GridSpec(2, 1,
                       #width_ratios=[1,1.2],
                       #height_ratios=[1,1]
                       )
# define vertical and horizontal spacing between panels
gs.update(wspace=0.3,hspace=0.4)

# possibly change outer margins of the figure
#plt.subplots_adjust(left=0.14, right=0.92, top=0.92, bottom=0.18)

gssub0 = gridspec.GridSpecFromSubplotSpec(1, len(nImgsBD) , subplot_spec=gs[0],wspace=0.05)

startN = 0
nInter = 0
for n in range(len(nImgsBD)):
    ax0 = plt.subplot(gssub0[n])
    for i in range(20): #range(len(intersectionROIs)):
        if intersectionROIs[i][8]>0.5:
            ax0.plot(timeStampsBD[n][1][:,1],FBD[intersectionROIs[i][0]][startN:(startN+len(timeStampsBD[n][1]))]+i*50,lw=0.2)
            if n==0:
                nInter+=1
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.spines['bottom'].set_position(('outward', 10))
    ax0.xaxis.set_ticks_position('bottom')
    if n==0:
        ax0.spines['left'].set_position(('outward', 10))
        ax0.yaxis.set_ticks_position('left')
        ax0.set_xlabel('time (s)')
    else:
        ax0.spines['left'].set_visible(False)
        ax0.yaxis.set_visible(False) #set_ticks_position('left')
    startN+=len(timeStampsBD[n][1])

print('intersection with more than 0.5 ',nInter)


gssub0 = gridspec.GridSpecFromSubplotSpec(1, len(nImgsAD) , subplot_spec=gs[1],wspace=0.05) #####################
#ax0.set_title('after drug')

startN = 0
nInter = 0
for n in range(len(nImgsAD)):
    ax0 = plt.subplot(gssub0[n])
    for i in range(20): #range(len(intersectionROIs)):
        if intersectionROIs[i][8]>0.5:
            ax0.plot(timeStampsAD[n][1][:,1],FAD[intersectionROIs[i][1]][startN:(startN+len(timeStampsAD[n][1]))]+i*50,lw=0.2)
            if n==0:
                nInter+=1
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.spines['bottom'].set_position(('outward', 10))
    ax0.xaxis.set_ticks_position('bottom')
    if n==0:
        ax0.spines['left'].set_position(('outward', 10))
        ax0.yaxis.set_ticks_position('left')
        ax0.set_xlabel('time (s)')
    else:
        ax0.spines['left'].set_visible(False)
        ax0.yaxis.set_visible(False) #set_ticks_position('left')
    startN+=len(timeStampsAD[n][1])


#plt.savefig(figOutDir+'FluorescenceTraces_%s.pdf' % animalID)
plt.show()