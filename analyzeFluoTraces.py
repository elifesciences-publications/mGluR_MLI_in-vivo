import numpy as np
import matplotlib.pyplot as plt
import pickle
import pdb
from ScanImageTiffReader import ScanImageTiffReader
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
from scipy import integrate
from scipy import stats

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
timeStampsBD = pickle.load(open( dataOutDir + 'timeStampsBeforeDrug_%s.p' % animalID, 'rb'))
timeStampsAD = pickle.load(open( dataOutDir + 'timeStampsAfterDrug_%s.p' % animalID, 'rb'))
#pdb.set_trace()

#################################################################
# calculate dF/F0 and integral

minOverlap = 0.5
startN  = 0
baselineTime = 6.
activityTimes = [10.,25.]

analysisBAndAD = []
nIntersections = 0
for i in range(len(intersectionROIs)):
    if intersectionROIs[i][8]>minOverlap:
        startNBD = 0
        beforeDrugAnalysis = []
        for n in range(len(timeStampsBD)):
            tBD = timeStampsBD[n][1][:,1]
            FBDsingle = FBD[intersectionROIs[i][0]][startNBD:(startNBD+len(timeStampsBD[n][1][:,1]))] - 0.7*FneuBD[intersectionROIs[i][0]][startNBD:(startNBD+len(timeStampsBD[n][1][:,1]))]
            baselineMask = tBD < baselineTime
            F0BD = np.mean(FBDsingle[baselineMask])
            deltaBD = (FBDsingle-F0BD)/F0BD
            activityMask = (tBD > activityTimes[0]) & (tBD < activityTimes[1])
            integralBD = integrate.simps(deltaBD, tBD)
            #integralBD = np.sum(deltaBD)
            FActivityBD = np.mean(FBDsingle[activityMask])
            beforeDrugAnalysis.append([n,integralBD,FActivityBD,F0BD,np.column_stack((tBD,FBDsingle,deltaBD))])
            #iE.append(integralBD)
            startNBD += len(timeStampsBD[n][1])
            if n==0:
                nIntersections+=1
        startNAD = 0
        afterDrugAnalysis = []
        for n in range(len(timeStampsAD)):
            tAD = timeStampsAD[n][1][:,1]
            FADsingle = FAD[intersectionROIs[i][1]][startNAD:(startNAD+len(timeStampsAD[n][1][:,1]))] - 0.7*FneuAD[intersectionROIs[i][1]][startNAD:(startNAD+len(timeStampsAD[n][1][:,1]))]
            baselineMask = tAD < baselineTime
            F0AD = np.mean(FADsingle[baselineMask])
            deltaAD = (FADsingle-F0AD)/F0AD
            activityMask = (tAD > activityTimes[0]) & (tAD < activityTimes[1])
            integralAD = integrate.simps(deltaAD, tAD)
            #integralAD = np.sum(deltaAD)
            FActivityAD = np.mean(FADsingle[activityMask])
            afterDrugAnalysis.append([n, integralAD, FActivityAD, F0AD, np.column_stack((tAD, FADsingle, deltaAD))])
            startNAD += len(timeStampsAD[n][1])
        analysisBAndAD.append([intersectionROIs[i],beforeDrugAnalysis,afterDrugAnalysis])

print('total number of intersection ROIs :', nIntersections)
pickle.dump(analysisBAndAD, open( dataOutDir + 'analysisBeforeAndAfterDrug_%s.p' % animalID, 'wb' ) )


##################################################################
# Show final results

fig_width = 20 # width in inches
fig_height = 15  # height in inches
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
# create figure instance
fig = plt.figure()

# define sub-panel grid and possibly width and height ratios
gs = gridspec.GridSpec(4, 1,
                       #width_ratios=[1,1.2],
                       #height_ratios=[1,1]
                       )
# define vertical and horizontal spacing between panels
gs.update(wspace=0.3,hspace=0.35)

# possibly change outer margins of the figure
plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)

gssub0 = gridspec.GridSpecFromSubplotSpec(1, len(timeStampsBD) , subplot_spec=gs[0],wspace=0.05)

# creating of figure instances
figsBD = []
for n in range(len(analysisBAndAD[0][1])):
   ax0 = plt.subplot(gssub0[n])
   ax0.spines['top'].set_visible(False)
   ax0.spines['right'].set_visible(False)
   ax0.spines['bottom'].set_position(('outward', 10))
   ax0.xaxis.set_ticks_position('bottom')
   if n == 0:
       ax0.spines['left'].set_position(('outward', 10))
       ax0.yaxis.set_ticks_position('left')
       ax0.set_xlabel('time (s)')
   else:
       ax0.spines['left'].set_visible(False)
       ax0.yaxis.set_visible(False)  # set_ticks_position('left')
   figsBD.append(ax0)

for i in range(10): #len(analysisBAndAD)):
    for n in range(len(analysisBAndAD[i][1])):
        #pdb.set_trace()
        figsBD[n].plot(analysisBAndAD[i][1][n][4][:,0],analysisBAndAD[i][1][n][4][:,1]+i*70,lw=0.2)


gssub1 = gridspec.GridSpecFromSubplotSpec(1, len(timeStampsAD) , subplot_spec=gs[1],wspace=0.05) #####################


# creating of figure instances
figsAD = []
for n in range(len(analysisBAndAD[0][2])):
   ax0 = plt.subplot(gssub1[n])
   ax0.spines['top'].set_visible(False)
   ax0.spines['right'].set_visible(False)
   ax0.spines['bottom'].set_position(('outward', 10))
   ax0.xaxis.set_ticks_position('bottom')
   if n == 0:
       ax0.spines['left'].set_position(('outward', 10))
       ax0.yaxis.set_ticks_position('left')
       ax0.set_xlabel('time (s)')
   else:
       ax0.spines['left'].set_visible(False)
       ax0.yaxis.set_visible(False)  # set_ticks_position('left')
   figsAD.append(ax0)

for i in range(10): #len(analysisBAndAD)):
    for n in range(len(analysisBAndAD[i][2])):
        #pdb.set_trace()
        figsAD[n].plot(analysisBAndAD[i][2][n][4][:,0],analysisBAndAD[i][2][n][4][:,1]+i*70,lw=0.2)


gssub2 = gridspec.GridSpecFromSubplotSpec(1, 4 , subplot_spec=gs[2],wspace=0.4) #####################

ax0 = plt.subplot(gssub2[0])
ax0.set_title('activity integral')
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_position(('outward', 10))
ax0.xaxis.set_ticks_position('bottom')
ax0.spines['left'].set_position(('outward', 10))
ax0.yaxis.set_ticks_position('left')
ax0.set_xlabel('recording number')

ax1 = plt.subplot(gssub2[1])
ax1.set_title('mean activity during locomotion')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.xaxis.set_ticks_position('bottom')
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.set_xlabel('recording number')

ax2 = plt.subplot(gssub2[2])
ax2.set_title('F0')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_position(('outward', 10))
ax2.xaxis.set_ticks_position('bottom')
ax2.spines['left'].set_position(('outward', 10))
ax2.yaxis.set_ticks_position('left')
ax2.set_xlabel('recording number')

trialsBeforeDrug = len(analysisBAndAD[0][1])
trialsAfterDrug = len(analysisBAndAD[0][2])

integralEv = []
activityEv = []
F0Ev = []

mColors = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9']
for i in range(len(analysisBAndAD)):
    #if analysisBAndAD[i][1][n][2] > 50:
    integral = []
    activity = []
    F0 = []
    for n in range(len(analysisBAndAD[i][1])):
        integral.append(analysisBAndAD[i][1][n][1])
        activity.append(analysisBAndAD[i][1][n][2])
        F0.append(analysisBAndAD[i][1][n][3])
    for m in range(len(analysisBAndAD[i][2])):
        integral.append(analysisBAndAD[i][2][m][1])
        activity.append(analysisBAndAD[i][2][m][2])
        F0.append(analysisBAndAD[i][2][m][3])

    integralEv.append([np.mean(integral[:trialsBeforeDrug]), np.mean(integral[trialsBeforeDrug:])])
    activityEv.append([np.mean(activity[:trialsBeforeDrug]), np.mean(activity[trialsBeforeDrug:])])
    F0Ev.append([np.mean(F0[:trialsBeforeDrug]), np.mean(F0[trialsBeforeDrug:])])

    if i<10:
        ax0.plot(range(1,trialsBeforeDrug+trialsAfterDrug+1),integral,'o-',ms=2,c=mColors[i%10],alpha=0.3)
        ax1.plot(range(1,trialsBeforeDrug+trialsAfterDrug+1),activity,'o-',ms=2,c=mColors[i%10],alpha=0.3)
        ax2.plot(range(1,trialsBeforeDrug+trialsAfterDrug+1),F0,'o-',ms=2,c=mColors[i%10],alpha=0.3)

        ax0.plot(range(trialsBeforeDrug + 1, trialsBeforeDrug + trialsAfterDrug + 1), integral[trialsBeforeDrug:], 'o-', ms=2,c=mColors[i%10],)
        ax1.plot(range(trialsBeforeDrug + 1, trialsBeforeDrug + trialsAfterDrug + 1), activity[trialsBeforeDrug:], 'o-', ms=2,c=mColors[i%10],)
        ax2.plot(range(trialsBeforeDrug + 1, trialsBeforeDrug + trialsAfterDrug + 1), F0[trialsBeforeDrug:], 'o-', ms=2,c=mColors[i%10],)



gssub3 = gridspec.GridSpecFromSubplotSpec(1, 4 , subplot_spec=gs[3],wspace=0.4) #####################

integralEv = np.asarray(integralEv)
activityEv = np.asarray(activityEv)
F0Ev = np.asarray(F0Ev)

ax0 = plt.subplot(gssub3[0])
a = [-10,60]
ax0.plot(a,a,c='0.5')
ax0.plot(integralEv[:,0],integralEv[:,1],'o',ms=3,alpha=0.5)

ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_position(('outward', 10))
ax0.xaxis.set_ticks_position('bottom')
ax0.spines['left'].set_position(('outward', 10))
ax0.yaxis.set_ticks_position('left')
ax0.set_xlabel('before drug')
ax0.set_ylabel('after drug')


ax1 = plt.subplot(gssub3[1])
a = [0,100]
ax1.plot(a,a,c='0.5')
ax1.plot(activityEv[:,0],activityEv[:,1],'o',ms=2,alpha=0.5)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.xaxis.set_ticks_position('bottom')
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.set_xlabel('before drug')
ax1.set_ylabel('after drug')

ax2 = plt.subplot(gssub3[2])
a = [0,100]
ax2.plot(a,a,c='0.5')
ax2.plot(F0Ev[:,0],F0Ev[:,1],'o',ms=2,alpha=0.5)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_position(('outward', 10))
ax2.xaxis.set_ticks_position('bottom')
ax2.spines['left'].set_position(('outward', 10))
ax2.yaxis.set_ticks_position('left')
ax2.set_xlabel('before drug')
ax2.set_ylabel('after drug')

print('paired t-test for integral : ',stats.ttest_rel(integralEv[:,0],integralEv[:,1]))
print('paired t-test for mean activity : ',stats.ttest_rel(activityEv[:,0],activityEv[:,1]))
print('paired t-test for F0 : ',stats.ttest_rel(F0Ev[:,0],F0Ev[:,1]))

plt.savefig(figOutDir+'FluorescenceTraces_%s.pdf' % animalID)
plt.show()