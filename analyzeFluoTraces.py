import numpy as np
import matplotlib.pyplot as plt
import scipy
import pickle
import pdb
from ScanImageTiffReader import ScanImageTiffReader
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
from scipy import integrate
from scipy import stats
import sys


from animalSettings import animalSettings
import graphicTools as GT

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

# moved to the parameter file
# baselineTime = 5.
# activityTimes = [10.,25.]

# waitTrialsAfterDrug = 3

mA = []
F0MeanBD = []
F0MeanAD = []
# determine least active period ###############################################################
meanActBD = []
meanActAD = []

for n in range(len(timeStampsBD)):  # loop over all recordings before drug delivery
    #
    twheelBD = wheelActivityBD[n][1][4]  # times of angles
    wActivityBD = (wheelActivityBD[n][1][3]) * 80. / 360.  # progress of animal angles conerted to  cm
    baselineWheelMask = twheelBD < aS.baselineTime
    meanAct = np.mean(np.fabs(wActivityBD[baselineWheelMask]))
    meanActBD.append(meanAct)

for n in range(len(timeStampsAD)):
    twheelAD = wheelActivityAD[n][1][4]
    wActivityAD = (wheelActivityAD[n][1][3]) * 80. / 360. # progress of animal angles conerted to  cm
    baselineWheelMask = twheelAD < aS.baselineTime
    meanAct = np.mean(np.fabs(wActivityAD[baselineWheelMask]))
    meanActAD.append(meanAct)

# determine with recording serves as F0
meanActBD = np.asarray(meanActBD)
meanActAD = np.asarray(meanActAD)

baselinePeriodBD = [] # np.zeros((2,len(timeStampsBD)))
baselinePeriodAD = [] #np.zeros((2,len(timeStampsAD)))

minTrialNumber = 1
# determine least active periods
for i in range(len(timeStampsBD)):
    for per in aS.baselinePeriodsBD :
        if per[0]<=i and i<=per[1]:
            threshold = 0.
            slowRecMask = [0,0] # first mask which should simply add to zero
            while sum(slowRecMask)<minTrialNumber:
                slowRecMask = (meanActBD < threshold) & ((np.arange(len(timeStampsBD)) >= per[0]) & (np.arange(len(timeStampsBD)) <= per[1]))  #meanActBD[per[0]:per[1]+1] < 100.
                threshold+=0.0001
            #print('slow BDslowRecMask)
            baselinePeriodBD.append([i,slowRecMask])

for i in range(len(timeStampsAD)):
    for per in aS.baselinePeriodsAD :
        if per[0]<=i and i<=per[1]:
            threshold = 0.
            slowRecMask = [0,0] # first mask which should simply add to zero
            while sum(slowRecMask)<minTrialNumber:
                slowRecMask = (meanActAD < threshold) & ((np.arange(len(timeStampsAD)) >= per[0]) & (np.arange(len(timeStampsAD)) <= per[1]))  #meanActBD[per[0]:per[1]+1] < 100.
                threshold+=0.0001
            #print(slowRecMask)
            baselinePeriodAD.append([i,slowRecMask])

print('BD :', baselinePeriodBD)
print('AD :', baselinePeriodAD)
#pdb.set_trace()

# determine F0 for non-active trials
nIntersection = 0
#F0allROIsBD = []
#F0allROIsAD = []
for i in range(len(intersectionROIs)): # loop over all ROIs
    if intersectionROIs[i][8]>aS.minOverlap: # only consider ROIs with minimal overlap
        startNBD = 0
        beforeDrugAnalysis = []
        F0perROIBD = []
        for n in range(len(timeStampsBD)): # loop over all recordings
            tBD = timeStampsBD[n][1][:, 1]
            FBDsingle = FBD[intersectionROIs[i][0]][startNBD:(startNBD + len(timeStampsBD[n][1][:, 1]))] - 0.7 * FneuBD[intersectionROIs[i][0]][startNBD:(startNBD + len(timeStampsBD[n][1][:, 1]))]
            baselineMask = tBD < aS.baselineTime  # only take period before motorization
            F0MBD = np.mean(FBDsingle[baselineMask]) # calculate mean of fluorescence during that period
            startNBD += len(timeStampsBD[n][1])
            F0perROIBD.append([n,F0MBD])
        F0perROIBD = np.asarray(F0perROIBD)
        #F0allROIsBD.append(F0perROIBD[:,1])
        F0perRecBD = []
        for n in range(len(timeStampsBD)): # loop over all recordings again to calculate mean F0 per recording
            #pdb.set_trace()
            F0m = np.mean(F0perROIBD[baselinePeriodBD[n][1]],axis=0)[1]
            F0perRecBD.append([n,F0m])
        #pdb.set_trace()
        F0MeanBD.append([i, nIntersection, F0perRecBD,F0perROIBD[:,1]])
        ######## after drug
        startNAD = 0
        F0perROIAD = []
        for n in range(len(timeStampsAD)):
            tAD = timeStampsAD[n][1][:, 1]
            FADsingle = FAD[intersectionROIs[i][1]][startNAD:(startNAD + len(timeStampsAD[n][1][:, 1]))] - 0.7 * FneuAD[intersectionROIs[i][1]][startNAD:(startNAD + len(timeStampsAD[n][1][:, 1]))]
            baselineMask = tAD < aS.baselineTime
            F0MAD = np.mean(FADsingle[baselineMask])
            startNAD += len(timeStampsAD[n][1])
            F0perROIAD.append([n,F0MAD])
        F0perROIAD = np.asarray(F0perROIAD)
        #F0allROIsAD.append(F0perROIAD[:, 1])
        F0perRecAD = []
        for n in range(len(timeStampsAD)):  # loop over all recordings again to calculate mean F0 per recording
            #pdb.set_trace()
            F0m = np.mean(F0perROIAD[baselinePeriodAD[n][1]],axis=0)[1]
            F0perRecAD.append([n, F0m])
        F0MeanAD.append([i, nIntersection, F0perRecAD,F0perROIAD[:, 1]])
        #####
        nIntersection+=1


#pdb.set_trace()
analysisBAndAD = []
nIntersection = 0
for i in range(len(intersectionROIs)):
    useInAnalysis = True
    if intersectionROIs[i][8]>aS.minOverlap:
        startNBD = 0
        beforeDrugAnalysis = []
        for n in range(len(timeStampsBD)):
            tBD = timeStampsBD[n][1][:,1]
            FBDsingle = FBD[intersectionROIs[i][0]][startNBD:(startNBD+len(timeStampsBD[n][1][:,1]))] - 0.7*FneuBD[intersectionROIs[i][0]][startNBD:(startNBD+len(timeStampsBD[n][1][:,1]))]
            #baselineMask = tBD < aS.baselineTime
            #F0BD = np.mean(FBDsingle[baselineMask])
            #F0BD = F0MeanBD[nIntersection]
            #if nIntersection != F0BD[1]:
            #    print('problem')
            #pdb.set_trace()
            F0Container =  F0MeanBD[nIntersection]
            if F0Container[1] == nIntersection:
                F0 = F0Container[2][n][1]
            else:
                print('problem')
            deltaBD = (FBDsingle - F0)
            deltaBDF0 = (FBDsingle-F0)/np.abs(F0)
            if F0<aS.F0minimum:
                print('excluded BD F0:',F0)
                useInAnalysis = False
            activityMask = (tBD > aS.activityTimes[0]) & (tBD < aS.activityTimes[1])
            integralBD = integrate.simps(deltaBDF0[activityMask], tBD[activityMask])
            #integralBD = np.sum(deltaBD)
            #FActivityBD = np.mean(deltaBDF0[activityMask])
            FActivityBD = np.percentile(deltaBDF0[activityMask],75)
            beforeDrugAnalysis.append([n,integralBD,FActivityBD,F0,np.column_stack((tBD,FBDsingle,deltaBDF0,deltaBD)),F0Container])
            #iE.append(integralBD)
            startNBD += len(timeStampsBD[n][1])
        startNAD = 0
        afterDrugAnalysis = []
        for n in range(len(timeStampsAD)):
            tAD = timeStampsAD[n][1][:,1]
            #pdb.set_trace()
            FADsingle = FAD[intersectionROIs[i][1]][startNAD:(startNAD+len(timeStampsAD[n][1][:,1]))] - 0.7*FneuAD[intersectionROIs[i][1]][startNAD:(startNAD+len(timeStampsAD[n][1][:,1]))]
            #baselineMask = tAD < aS.baselineTime
            #F0AD = F0MeanAD[nIntersections]
            #if nIntersections != F0AD[1]:
            #    print('problem')
            #F0AD = np.mean(FADsingle[baselineMask])
            F0Container = F0MeanAD[nIntersection]
            if F0Container[1] == nIntersection:
                F0 = F0Container[2][n][1]
            else:
                print('problem')
            deltaAD = (FADsingle-F0)
            deltaADF0 = (FADsingle-F0)/np.abs(F0)
            if F0<aS.F0minimum:
                print('excluded AD F0:',F0)
                useInAnalysis = False
            activityMask = (tAD > aS.activityTimes[0]) & (tAD < aS.activityTimes[1])
            #pdb.set_trace()
            #print(n)
            integralAD = integrate.simps(deltaADF0[activityMask],tAD[activityMask])
            #integralAD = np.sum(deltaAD[activityMask])
            #FActivityAD = np.mean(deltaADF0[activityMask])
            FActivityAD = np.percentile(deltaADF0[activityMask], 75)
            afterDrugAnalysis.append([n, integralAD, FActivityAD, F0, np.column_stack((tAD, FADsingle, deltaADF0, deltaAD)),F0Container])
            startNAD += len(timeStampsAD[n][1])
        if useInAnalysis:
            analysisBAndAD.append([intersectionROIs[i],beforeDrugAnalysis,afterDrugAnalysis])
        nIntersection += 1

#pdb.set_trace()
print('total number of intersection ROIs :', nIntersection)
pickle.dump(analysisBAndAD, open( dataOutDir + 'analysisBeforeAndAfterDrug_%s.p' % aS.animalID, 'wb' ) )

##################################################################
#GT.generateSingleTraceFigure(aS.animalID,timeStampsBD,timeStampsAD,analysisBAndAD,aS.baselinePeriodsBD,aS.baselinePeriodsAD,columnIdx=1)
#GT.generateSingleTraceFigure(aS.animalID,timeStampsBD,timeStampsAD,analysisBAndAD,aS.baselinePeriodsBD,aS.baselinePeriodsAD,columnIdx=2)
#GT.generateSingleTraceFigure(aS.animalID,timeStampsBD,timeStampsAD,analysisBAndAD,aS.baselinePeriodsBD,aS.baselinePeriodsAD,columnIdx=3)
#pdb.set_trace()
##################################################################
# Show final results

fig_width = 20 # width in inches
fig_height = 20  # height in inches
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
                       height_ratios=[1,1,2.5,1]
                       )
# define vertical and horizontal spacing between panels
gs.update(wspace=0.3,hspace=0.25)

# possibly change outer margins of the figure
plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)
plt.figtext(0.06, 0.97, '%s :   %s recs. w/o drug; %s recordings with drug' % (aS.animalID, len(analysisBAndAD[0][1]), len(analysisBAndAD[0][2])), clip_on=False, color='black', size=14)

gssub0 = gridspec.GridSpecFromSubplotSpec(2, len(timeStampsBD) , subplot_spec=gs[0],wspace=0.05,hspace=0.1,height_ratios=[0.2,1])

# creating of figure instances
figsBD = []
for n in range(len(analysisBAndAD[0][1])):
    ax0 = plt.subplot(gssub0[n])
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.spines['bottom'].set_visible(False)
    ax0.xaxis.set_visible(False)
    if n == 0:
        ax0.spines['left'].set_position(('outward', 10))
        ax0.yaxis.set_ticks_position('left')
        ax0.set_xlabel('time (s)')
    else:
        ax0.spines['left'].set_visible(False)
        ax0.yaxis.set_visible(False)  # set_ticks_position('left')


    ax1 = plt.subplot(gssub0[n+len(timeStampsBD)])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_position(('outward', 10))
    ax1.xaxis.set_ticks_position('bottom')
    if n == 0:
        ax1.spines['left'].set_position(('outward', 10))
        ax1.yaxis.set_ticks_position('left')
        ax1.set_xlabel('time (s)')
    else:
        ax1.spines['left'].set_visible(False)
        ax1.yaxis.set_visible(False)  # set_ticks_position('left')
    figsBD.append([ax0,ax1])

for i in range(10): #len(analysisBAndAD)):
    for n in range(len(analysisBAndAD[i][1])):
        #pdb.set_trace()
        if i==0:
            figsBD[n][0].plot(wheelActivityBD[n][1][4],wheelActivityBD[n][1][3])
        figsBD[n][1].plot(analysisBAndAD[i][1][n][4][:,0],analysisBAndAD[i][1][n][4][:,2]+i*2,lw=0.2)


gssub1 = gridspec.GridSpecFromSubplotSpec(2, len(timeStampsAD) , subplot_spec=gs[1],wspace=0.05,hspace=0.1,height_ratios=[0.2,1]) #####################


# creating of figure instances
figsAD = []
for n in range(len(analysisBAndAD[0][2])):
    ax0 = plt.subplot(gssub1[n])
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.spines['bottom'].set_visible(False)
    ax0.xaxis.set_visible(False)
    if n == 0:
        ax0.spines['left'].set_position(('outward', 10))
        ax0.yaxis.set_ticks_position('left')
        ax0.set_xlabel('time (s)')
    else:
        ax0.spines['left'].set_visible(False)
        ax0.yaxis.set_visible(False)  # set_ticks_position('left')

    ax1 = plt.subplot(gssub1[n+len(timeStampsAD)])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_position(('outward', 10))
    ax1.xaxis.set_ticks_position('bottom')
    if n == 0:
        ax1.spines['left'].set_position(('outward', 10))
        ax1.yaxis.set_ticks_position('left')
        ax1.set_xlabel('time (s)')
    else:
        ax1.spines['left'].set_visible(False)
        ax1.yaxis.set_visible(False)  # set_ticks_position('left')
    figsAD.append([ax0,ax1])

for i in range(10): #len(analysisBAndAD)):
    for n in range(len(analysisBAndAD[i][2])):
        #pdb.set_trace()
        if i==0:
            figsAD[n][0].plot(wheelActivityAD[n][1][4],wheelActivityAD[n][1][3])
        figsAD[n][1].plot(analysisBAndAD[i][2][n][4][:,0],analysisBAndAD[i][2][n][4][:,2]+i*2,lw=0.2)


gssub2 = gridspec.GridSpecFromSubplotSpec(2, 4 , subplot_spec=gs[2],wspace=0.4,hspace=0.4) #####################

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

ax3 = plt.subplot(gssub2[4])
ax3.set_title('activity integral \n trend corrected')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['bottom'].set_position(('outward', 10))
ax3.xaxis.set_ticks_position('bottom')
ax3.spines['left'].set_position(('outward', 10))
ax3.yaxis.set_ticks_position('left')
ax3.set_xlabel('recording number')

ax4 = plt.subplot(gssub2[5])
ax4.set_title('mean activity during locomotion \n trend corrected')
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.spines['bottom'].set_position(('outward', 10))
ax4.xaxis.set_ticks_position('bottom')
ax4.spines['left'].set_position(('outward', 10))
ax4.yaxis.set_ticks_position('left')
ax4.set_xlabel('recording number')


trialsBeforeDrug = len(analysisBAndAD[0][1])
trialsAfterDrug = len(analysisBAndAD[0][2])

integralEv = []
activityEv = []
F0Ev = []

mColors = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9']
totRecordings = len(analysisBAndAD[0][1]) + len(analysisBAndAD[0][2])
print('total recordings',totRecordings)
integralMean = np.zeros(totRecordings)
activityMean = np.zeros(totRecordings)

nChangeUP = 0
nChangeDOWN = 0
nNoChange = 0

nChangeUPAct = 0
nChangeDOWNAct = 0
nNoChangeAct = 0
if aS.animalID == 'animal#1_2':
    fit = 'linear'
else:
    fit = 'exponential'
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

    integral = np.asarray(integral)
    activity = np.asarray(activity)

    if fit=='linear':
        polycoeffsInt = scipy.polyfit(range(trialsBeforeDrug), integral[:trialsBeforeDrug], 1)
        integralFit = scipy.polyval(polycoeffsInt, range(1,trialsBeforeDrug+trialsAfterDrug+1))

        polycoeffsAct = scipy.polyfit(range(1, trialsBeforeDrug + 1), activity[:trialsBeforeDrug], 1)
        activityFit = scipy.polyval(polycoeffsAct, range(1, trialsBeforeDrug + trialsAfterDrug + 1))
    elif fit=='exponential':
        fitfunc = lambda p, x: p[0] + p[1] * np.exp(-x/p[2])
        errfunc = lambda p, x, y: fitfunc(p, x) - y + np.where(1000.,p[2]<0,0.)
        # fit a gaussian to the correlation function
        p0Int = [0.7*integral[0],0.3*integral[0],3]
        p1Int, success = scipy.optimize.leastsq(errfunc, p0Int, args=(np.arange(trialsBeforeDrug), integral[:trialsBeforeDrug]))
        integralFit = fitfunc(p1Int, np.arange(trialsBeforeDrug+trialsAfterDrug))
        p0Act = [0.7*activity[0],0.3*activity[0],3]
        p1Act, success = scipy.optimize.leastsq(errfunc, p0Act, args=(np.arange(trialsBeforeDrug), activity[:trialsBeforeDrug]))
        activityFit = fitfunc(p1Act, np.arange(trialsBeforeDrug+trialsAfterDrug))

    integralTrendCorrected = (integral - integralFit)#+np.mean(integral[:trialsBeforeDrug])
    activityTrendCorrected = (activity - activityFit)#+np.mean(activity[:trialsBeforeDrug])
    activityTrendCorrectedRenormalized = (activity - activityFit) +np.mean(activity[:trialsBeforeDrug])
    downInt = -1
    tttInt = stats.ttest_ind(integralTrendCorrected[:trialsBeforeDrug], integralTrendCorrected[(trialsBeforeDrug + aS.waitTrialsAfterDrug):])
    tttAct = stats.ttest_ind(activityTrendCorrected[:trialsBeforeDrug], activityTrendCorrected[(trialsBeforeDrug + aS.waitTrialsAfterDrug):])


    if tttInt[1] < 0.01 and (np.mean(integralTrendCorrected[:trialsBeforeDrug]) > np.mean(integralTrendCorrected[(trialsBeforeDrug + aS.waitTrialsAfterDrug):])):
        nChangeDOWN += 1
        downInt = 1
    elif tttInt[1] < 0.01 and (np.mean(integralTrendCorrected[:trialsBeforeDrug]) < np.mean(integralTrendCorrected[(trialsBeforeDrug + aS.waitTrialsAfterDrug):])):
        nChangeUP += 1
        downInt = 0
    else:
        nNoChange += 1
        downInt = 2

    downAct = -1
    if tttAct[1] < 0.01 and (np.mean(activityTrendCorrected[:trialsBeforeDrug]) > np.mean(activityTrendCorrected[(trialsBeforeDrug + aS.waitTrialsAfterDrug):])):
        nChangeDOWNAct += 1
        downAct = 1
    elif tttInt[1] < 0.01 and (np.mean(activityTrendCorrected[:trialsBeforeDrug]) < np.mean(activityTrendCorrected[(trialsBeforeDrug + aS.waitTrialsAfterDrug):])):
        nChangeUPAct += 1
        downAct = 0
    else:
        nNoChangeAct += 1
        downAct = 2

    integralEv.append([np.mean(integral[:trialsBeforeDrug]), np.mean(integralTrendCorrected[(trialsBeforeDrug+aS.waitTrialsAfterDrug):]),tttInt[1],downInt])
    activityEv.append([np.mean(activity[:trialsBeforeDrug]), np.mean(activityTrendCorrected[(trialsBeforeDrug+aS.waitTrialsAfterDrug):]),tttAct[1],downAct,np.mean(activityTrendCorrectedRenormalized[:trialsBeforeDrug]),np.mean(activityTrendCorrectedRenormalized[(trialsBeforeDrug+aS.waitTrialsAfterDrug):])])
    F0Ev.append([np.mean(F0[:trialsBeforeDrug]), np.mean(F0[trialsBeforeDrug:])])

    if i==0:
        ax0.axhline(y=0,ls='--',c='0.5')
        ax1.axhline(y=0, ls='--', c='0.5')

    if i<10:
        ax0.plot(range(1,trialsBeforeDrug+trialsAfterDrug+1),integral,'o-',ms=2,c=mColors[i%10],alpha=0.3)
        ax1.plot(range(1,trialsBeforeDrug+trialsAfterDrug+1),activity,'o-',ms=2,c=mColors[i%10],alpha=0.3)
        ax2.plot(range(1,trialsBeforeDrug+trialsAfterDrug+1),F0,'o-',ms=2,c=mColors[i%10],alpha=0.3)
        #pdb.set_trace()
        ax2.plot(range(1,trialsBeforeDrug+1),analysisBAndAD[i][1][n][5][3],'D-',ms=2,c=mColors[i%10],alpha=0.3)


        ax0.plot(range(trialsBeforeDrug + 1, trialsBeforeDrug + trialsAfterDrug + 1), integral[trialsBeforeDrug:], 'o-', ms=2,c=mColors[i%10],)
        ax1.plot(range(trialsBeforeDrug + 1, trialsBeforeDrug + trialsAfterDrug + 1), activity[trialsBeforeDrug:], 'o-', ms=2,c=mColors[i%10],)
        ax2.plot(range(trialsBeforeDrug + 1, trialsBeforeDrug + trialsAfterDrug + 1), F0[trialsBeforeDrug:], 'o-', ms=2,c=mColors[i%10],)
        ax2.plot(range(trialsBeforeDrug + 1, trialsBeforeDrug + trialsAfterDrug + 1), analysisBAndAD[i][2][n][5][3], 'D-', ms=2, c=mColors[i % 10], alpha=0.3)

    if i==0:
        ax3.axhline(y=0,ls='--',c='0.7')
        ax4.axhline(y=0, ls='--', c='0.7')

    if i<10:
        ax3.plot(range(1, trialsBeforeDrug + trialsAfterDrug + 1), integral - integralFit, 'o-', ms=2, c=mColors[i % 10], alpha=0.3)
        ax4.plot(range(1, trialsBeforeDrug + trialsAfterDrug + 1), activity - activityFit, 'o-', ms=2, c=mColors[i % 10], alpha=0.3)

        ax3.plot(range(trialsBeforeDrug + 1, trialsBeforeDrug + trialsAfterDrug + 1), integral[trialsBeforeDrug:] - integralFit[trialsBeforeDrug:], 'o-', ms=2, c=mColors[i % 10], )
        ax4.plot(range(trialsBeforeDrug + 1, trialsBeforeDrug + trialsAfterDrug + 1), activity[trialsBeforeDrug:] - activityFit[trialsBeforeDrug:], 'o-', ms=2, c=mColors[i % 10], )



    integralMean += np.asarray(integral)/len(analysisBAndAD)
    activityMean += np.asarray(activity)/len(analysisBAndAD)



print('Integral of Rois changed down : ',nChangeDOWN)
print('Integral of Rois changed up : ',nChangeUP)
print('Integral Rois without change ',nNoChange)

print('Activity of Rois changed down : ',nChangeDOWNAct)
print('Activity of Rois changed up : ',nChangeUPAct)
print('Activity Rois without change ',nNoChangeAct)
ax0.plot(range(1,trialsBeforeDrug+trialsAfterDrug+1),integralMean,'o-',ms=2,c='black')
ax1.plot(range(1,trialsBeforeDrug+trialsAfterDrug+1),activityMean,'o-',ms=2,c='black')

integralEv = np.asarray(integralEv)
activityEv = np.asarray(activityEv)
F0Ev = np.asarray(F0Ev)

#ax0.plot(range(1,trialsBeforeDrug+trialsAfterDrug+1),np.mean(integralEv[:,0]),'o-',ms=2,c='black')
#ax1.plot(range(1,trialsBeforeDrug+trialsAfterDrug+1),np.mean(activityEv,axis=1),'o-',ms=2,c='black')





gssub3 = gridspec.GridSpecFromSubplotSpec(1, 4 , subplot_spec=gs[3],wspace=0.4) #####################



ax0 = plt.subplot(gssub3[0])
#a = [-10,60]
#ax0.plot(a,a,c='0.5')
ax0.axhline(y=0,ls='--',c='0.5')
for i in range(len(integralEv)):
    if integralEv[i,2]<0.01 and integralEv[i,3] == 1:
        ax0.plot(integralEv[i,0],integralEv[i,1],'o',ms=3,c='C0')
    elif integralEv[i,2]<0.01 and integralEv[i,3] == 0:
        ax0.plot(integralEv[i,0],integralEv[i,1],'o',ms=3,c='C1')
    else:
        ax0.plot(integralEv[i, 0], integralEv[i, 1], 'o',c='C2', ms=3,alpha=0.5)

#from scipy.stats import linregress
#mask = np.abs(integralEv[:,1])<100.
#out1 = linregress(integralEv[:,0][mask],integralEv[:,1][mask])
#print(out1)
#yPred = out1[1] + out1[0]*integralEv[:,0]
#ax0.plot(integralEv[:,0],yPred)

ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_position(('outward', 10))
ax0.xaxis.set_ticks_position('bottom')
ax0.spines['left'].set_position(('outward', 10))
ax0.yaxis.set_ticks_position('left')
ax0.set_xlabel('before drug')
ax0.set_ylabel('after drug')

if aS.intLims is not None:
    ax0.set_ylim(aS.intLims[0],aS.intLims[1])
    ax0.set_xlim(aS.intLims[2],aS.intLims[3])
#ax0.set_ylim(-60,40)

text11 = ax0.annotate('%s ROIs changed down \n%s ROIs changed up \n%s ROIs didn\'t change' %(nChangeDOWN,nChangeUP,nNoChange) , xy=(50,20), annotation_clip=False,
           xytext=None, textcoords='data',fontsize=10,
           arrowprops=None
           )

print('Integral of Rois changed down : ',nChangeDOWN)
print('Integral of Rois changed up : ',nChangeUP)
print('Integral Rois without change ',nNoChange)
ax1 = plt.subplot(gssub3[1])
#a = [0,3]
#ax1.plot(a,a,c='0.5')
ax1.axhline(y=0,ls='--',c='0.5')
for i in range(len(activityEv)):
    if activityEv[i,2]<0.01 and activityEv[i,3] == 1:
        ax1.plot(activityEv[i,0],activityEv[i,1],'o',ms=3,c='C0')
    elif activityEv[i,2]<0.01 and activityEv[i,3] == 0:
        ax1.plot(activityEv[i,0],activityEv[i,1],'o',ms=3,c='C1')
    else:
        ax1.plot(activityEv[i,0],activityEv[i,1],'o',c='C2',ms=2,alpha=0.5)

#ax1.plot(activityEv[:,0],activityEv[:,1],'o',ms=2,alpha=0.5)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.xaxis.set_ticks_position('bottom')
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.set_xlabel('mean before drug')
ax1.set_ylabel('change in activity after drug')
#ax1.set_ylim(-4,)

text11 = ax1.annotate('%s ROIs changed down \n %s ROIs changed up \n %s ROIs didn\'t change' %(nChangeDOWNAct,nChangeUPAct,nNoChangeAct) , xy=(3,2), annotation_clip=False,
           xytext=None, textcoords='data',fontsize=10,
           arrowprops=None
           )

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
ax2.set_xlabel('mean before drug')
ax2.set_ylabel('change after drug')

#print('paired t-test for integral : ',stats.ttest_rel(integralEv[:,0],integralEv[:,1]))
#print('paired t-test for mean activity : ',stats.ttest_rel(activityEv[:,0],activityEv[:,1]))
#print('paired t-test for F0 : ',stats.ttest_rel(F0Ev[:,0],F0Ev[:,1]))

intMask = integralEv[:,0]>10
actMask = activityEv[:,0]>0.7

print('t-test for integral : ',stats.ttest_1samp(integralEv[:,1][intMask],0.))
print('t-test for mean activity : ',stats.ttest_1samp(activityEv[:,1][actMask],0.))
print('paired t-test for F0 : ',stats.ttest_rel(F0Ev[:,0],F0Ev[:,1]))

plt.savefig(figOutDir+'FluorescenceTraces_%s.pdf' % aS.animalID)
#plt.savefig(figOutDir+'FluorescenceTraces_%s.png' % animalID)
#plt.show()

allData = [trialsBeforeDrug,trialsAfterDrug,integral,activity,integralEv,activityEv,F0Ev]
pickle.dump(allData, open( dataOutDir + 'activityChanges_%s.p' % aS.animalID, 'wb' ) )
