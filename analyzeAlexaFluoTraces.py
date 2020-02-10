import numpy as np
import pickle
import sys
import re
import datetime
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
import pdb

from ScanImageTiffReader import ScanImageTiffReader


from animalSettings import animalSettings
import graphicTools as GT

def readTimeStampOfRecording(tiffFile,nFrame):
   desc = ScanImageTiffReader(tiffFile).description(nFrame)
   keyWordIdx = desc.find('epoch')
   dateString = re.split('\[|\]', desc[keyWordIdx:])
   dateIdv = dateString[1].split()
   #print('dateString',dateString)
   #print('dateIdv',dateIdv)
   unixStartTime = int(datetime.datetime(int(dateIdv[0]),int(dateIdv[1]),int(dateIdv[2]),int(dateIdv[3]),int(dateIdv[4]),int(float(dateIdv[5]))).strftime('%s'))
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
   local_time = datetime.datetime.fromtimestamp(unixFrameTime)
   #print(type(dt))
   print(local_time.strftime("%Y-%m-%d %H:%M:%S.%f%z"))
   DDate = local_time.strftime("%Y-%m-%d")
   TTime = local_time.strftime("%H:%M:%S.%f%z")
   #return ([frameNumberAcquisition,frameTimestamps-acqTriggerTimestamps,acquisitionNumbers,unixStartTime,unixFrameTime,frameTimestamps,acqTriggerTimestamps])
   return ([frameNumberAcquisition,int(local_time.strftime("%Y")),int(local_time.strftime("%m")),int(local_time.strftime("%d")),int(local_time.strftime("%H")),int(local_time.strftime("%M")),float(local_time.strftime("%S.%f%z"))])



try:
    anim = sys.argv[1]
except IndexError:
    anim = 'animal#3'
else:
    pass

aS = animalSettings(anim)
aS.loadSettings()

print('animal : ',aS.animalID)


#########################################################
# set parameters

dataOutDir = 'dataOutput/'
figOutDir = 'figureOutput/'

# animalID = 'animal#1'
# #baseDir = '/media/labo_rw/JINmGluR/'
# baseDir = '/home/mgraupe/'
# beforeDrugDir = '2019.08.15_002_animal#1/animal#1_00000-00013/'
# afterDrugDir = '2019.08.15_002_animal#1/animal#1_00015-20_23-28/'
# suite2pDir = 'suite2p/plane0/'

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

FBD = np.load(aS.baseDir + aS.beforeDrugDir + aS.suite2pDir + 'F.npy')
FneuBD = np.load(aS.baseDir + aS.beforeDrugDir + aS.suite2pDir + 'Fneu.npy')
FAD = np.load(aS.baseDir + aS.afterDrugDir + aS.suite2pDir + 'F.npy')
FneuAD = np.load(aS.baseDir + aS.afterDrugDir + aS.suite2pDir + 'Fneu.npy')

#################################################################
# read ROI matching results
intersectionROIs = pickle.load(open( dataOutDir + 'ROIintersections_%s.p' % aS.animalID, 'rb'))

# read time stamps
timeStampsBD = pickle.load(open( dataOutDir + 'timeStampsBeforeDrug_%s.p' % aS.animalID, 'rb'))
timeStampsAD = pickle.load(open( dataOutDir + 'timeStampsAfterDrug_%s.p' % aS.animalID, 'rb'))
#pdb.set_trace()
# read wheel activity
wheelActivityBD = pickle.load(open( dataOutDir + 'wheelActivityBeforeDrug_%s.p' % aS.animalID, 'rb'))
wheelActivityAD = pickle.load(open( dataOutDir + 'wheelActivityAfterDrug_%s.p' % aS.animalID, 'rb'))

#pdb.set_trace()
#################################################################
# calculate dF/F0 and integral


#linkDir = '/media/labo_rw/JINmGluR/2019.08.15_002_animal#1/animal#1_00036-50/rawImages' #/home/mgraupe/experiments/in_vivo_cerebellum_walking/LocoRungs/testScripts/suite2pOutput'
#baseDir = '/media/invivodata2/altair_data/dataMichael'
#dataFolder = '2019.08.15_002'
#fileBaseName = 'animal#1_'
#startN = 36
#endN   = 50
#imgList = None # [1,2,3,4,5,6,8,9,10,11,12,13]
#startTime = [17,6,0]

if aS.animalID == 'animal#1_2':
    aIDFile = 'animal#1'
else:
    aIDFile = aS.animalID

FAlexaData = []
for i in aS.alexaImgList:
    # copy file to rawImages Direcotry
    print(i)
    tiffStack = '%s/%s/%s_%05d.tif' % (aS.rawDataDir,aS.dataFolder,aIDFile,i)
    #myCmd = 'cp %s/%s/%s%05d.tif %s/' % (baseDir,dataFolder,fileBaseName,i,linkDir)
    #os.system(myCmd)
    # extract time-stamps
    #timeStamps = []
    #f = '%s/%s%05d.tif' % (linkDir,fileBaseName,i)
    FData = ScanImageTiffReader(tiffStack).data()
    fN = np.shape(FData)[0]
    ttt = []
    for n in range(fN):
        timeStamp = readTimeStampOfRecording(tiffStack,n)
        timeDiff = (timeStamp[4] - aS.alexaStartTime[0])*60. + (timeStamp[5] - aS.alexaStartTime[1]) + (timeStamp[6] - aS.alexaStartTime[2])/60.
        #FAlexaData.append([readTimeStampOfRecording(tiffStack,n),np.average(FData[n])])
        ttt.append(timeDiff)
    FAlexaData.append([i,np.average(timeDiff),np.average(FData),np.std(np.average(FData,axis=(1,2)))])
    #pdb.set_trace()
    #pdb.set_trace()
    #timeStamps = np.asarray(timeStamps)
    #igorFileName = f[:-4]+'_timeStamps.csv'
    #np.savetxt(igorFileName,timeStamps,delimiter=',')
    #sio.savemat(matlabFileName, mdict={'timeStamps': timeStamps})

FAlexaData = np.asarray(FAlexaData)
# moved to the parameter file
# baselineTime = 5.
# activityTimes = [10.,25.]

# waitTrialsAfterDrug = 3

#pdb.set_trace()
#print('total number of intersection ROIs :', nIntersection)
pickle.dump(FAlexaData, open( dataOutDir + 'analysisOfAlexaRecordings_%s.p' % aS.animalID, 'wb' ) )

##################################################################
#GT.generateSingleTraceFigure(aS.animalID,timeStampsBD,timeStampsAD,analysisBAndAD,aS.baselinePeriodsBD,aS.baselinePeriodsAD,columnIdx=1)
#GT.generateSingleTraceFigure(aS.animalID,timeStampsBD,timeStampsAD,analysisBAndAD,aS.baselinePeriodsBD,aS.baselinePeriodsAD,columnIdx=2)
#GT.generateSingleTraceFigure(aS.animalID,timeStampsBD,timeStampsAD,analysisBAndAD,aS.baselinePeriodsBD,aS.baselinePeriodsAD,columnIdx=3)
#pdb.set_trace()
##################################################################
# Show final results

fig_width = 8 # width in inches
fig_height = 7  # height in inches
fig_size =  [fig_width,fig_height]
params = {'axes.labelsize': 12,
          'axes.titlesize': 12,
          'font.size': 12,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
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
gs = gridspec.GridSpec(1, 1
                       #width_ratios=[1,1.2],
                       #height_ratios=[1,1,2.5,1]
                       )
# define vertical and horizontal spacing between panels
gs.update(wspace=0.3,hspace=0.25)

# possibly change outer margins of the figure
plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)
plt.figtext(0.06, 0.97, '%s :   %s recordings of each %s images with Alex 594' % (aS.animalID, len(aS.alexaImgList), fN), clip_on=False, color='black', size=14)

#gssub0 = gridspec.GridSpecFromSubplotSpec(2, len(timeStampsBD) , subplot_spec=gs[0],wspace=0.05,hspace=0.1,height_ratios=[0.2,1])
ax0 = plt.subplot(gs[0])
# creating of figure instances
ax0.errorbar(FAlexaData[:,1],FAlexaData[:,2],FAlexaData[:,3],fmt='o-')
#ax0.set_title('Alexa 594 fluorescence')
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_position(('outward', 10))
ax0.xaxis.set_ticks_position('bottom')
ax0.spines['left'].set_position(('outward', 10))
ax0.yaxis.set_ticks_position('left')

#ax0.set_ylim(0,5)

ax0.set_xlabel('time since Alexa 594 application (min)')
ax0.set_ylabel('average fluorescence per image')


#plt.show()
plt.savefig(figOutDir+'Alexa594FluorescenceTraces_%s.pdf' % aS.animalID)
#plt.savefig(figOutDir+'FluorescenceTraces_%s.png' % animalID)
plt.show()

#allData = [trialsBeforeDrug,trialsAfterDrug,integral,activity,integralEv,activityEv,F0Ev]
#pickle.dump(allData, open( dataOutDir + 'activityChanges_%s.p' % aS.animalID, 'wb' ) )
