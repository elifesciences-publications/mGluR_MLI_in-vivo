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
from scipy.stats import linregress

from animalSettings import animalSettings
import graphicTools as GT

def movingAverage(x, w):
    return np.convolve(x, np.ones(w), 'same')/w

def layoutOfPanel(ax,xLabel=None,yLabel=None,Leg=None,xyInvisible=[False,False]):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #
    if xyInvisible[0]:
        ax.spines['bottom'].set_visible(False)
        ax.xaxis.set_visible(False)
    else:
        ax.spines['bottom'].set_position(('outward', 10))
        ax.xaxis.set_ticks_position('bottom')
    #
    if xyInvisible[1]:
        ax.spines['left'].set_visible(False)
        ax.yaxis.set_visible(False)
    else:
        ax.spines['left'].set_position(('outward', 10))
        ax.yaxis.set_ticks_position('left')


    if xLabel != None :
        ax.set_xlabel(xLabel,fontsize=20)

    if yLabel != None :
        ax.set_ylabel(yLabel,fontsize=20)

    if Leg != None :
        ax.legend(loc=Leg[0], frameon=False)
        if len(Leg)>1 :
            legend = ax.get_legend()  # plt.gca().get_legend()
            ltext = legend.get_texts()
            plt.setp(ltext, fontsize=Leg[1])

################################################


animal = 'animal#2'

aS = animalSettings(animal)
aS.loadSettings()


print('example based on animal : ',aS.animalID)

#########################################################
# set parameters

dataOutDir = 'dataOutput/'
figOutDir = 'publicationFigures/'

##########################################################
# read data and determine principal parameters
# statBD = np.load(aS.baseDir + aS.beforeDrugDir + aS.suite2pDir + 'stat.npy')
opsBD = np.load(aS.baseDir + aS.beforeDrugDir + aS.suite2pDir + 'ops.npy').item()
# iscellBD = np.load(aS.baseDir + aS.beforeDrugDir + aS.suite2pDir +'iscell.npy')
# ncellsBD = len(iscellBD)
#
# statAD = np.load(aS.baseDir + aS.afterDrugDir + aS.suite2pDir + 'stat.npy')
# opsAD = np.load(aS.baseDir + aS.afterDrugDir + aS.suite2pDir + 'ops.npy').item()
# iscellAD = np.load(aS.baseDir + aS.afterDrugDir + aS.suite2pDir +'iscell.npy')
# ncellsAD = len(iscellAD)
#
# FBD = np.load(aS.baseDir + aS.beforeDrugDir + aS.suite2pDir + 'F.npy')
# FneuBD = np.load(aS.baseDir + aS.beforeDrugDir + aS.suite2pDir + 'Fneu.npy')
# FAD = np.load(aS.baseDir + aS.afterDrugDir + aS.suite2pDir + 'F.npy')
# FneuAD = np.load(aS.baseDir + aS.afterDrugDir + aS.suite2pDir + 'Fneu.npy')

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

##################################################################
# read data for figure

#np.save('imageBD.npy',opsBD['meanImgE'])
imgBD = np.load('imageBD.npy')
print('saved')
imBD =  imgBD[aS.cutOffXPub:-aS.cutOffXPub,aS.cutOffYPub:-(aS.cutOffYPub)]
print('image size : ', np.shape(imBD))
analysisBAndAD = pickle.load(open( dataOutDir + 'analysisBeforeAndAfterDrug_%s.p' % aS.animalID, 'rb' ) )

ROIIDs = [0,1,31,3,17,5,19]
recordingsBD = [0,3,5]
recordingsAD  = [0,4,7]

trialsBeforeDrug = len(analysisBAndAD[0][1])
trialsAfterDrug = len(analysisBAndAD[0][2])

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
params = {'axes.labelsize': 18,
          'axes.titlesize': 16,
          'font.size': 16,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
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
                       height_ratios=[2,1]
                       )
# define vertical and horizontal spacing between panels
gs.update(wspace=0.3,hspace=0.17)

# possibly change outer margins of the figure
plt.subplots_adjust(left=0.07, right=0.97, top=0.97, bottom=0.07)
#plt.figtext(0.06, 0.97, '%s :   %s recs. w/o drug; %s recordings with drug' % (aS.animalID, len(analysisBAndAD[0][1]), len(analysisBAndAD[0][2])), clip_on=False, color='black', size=14)
plt.figtext(0.03, 0.975, 'A',clip_on=False,color='black', weight='bold',size=28)

plt.figtext(0.38, 0.975, 'C',clip_on=False,color='black', weight='bold',size=28)

plt.figtext(0.03, 0.68, 'B',clip_on=False,color='black', weight='bold',size=28)
plt.figtext(0.38, 0.68, 'D',clip_on=False,color='black', weight='bold',size=28)
#plt.figtext(0.675, 0.616, 'E',clip_on=False,color='black', weight='bold',size=28)

plt.figtext(0.03, 0.36, 'E',clip_on=False,color='black', weight='bold',size=28)
plt.figtext(0.32, 0.36, 'F',clip_on=False,color='black', weight='bold',size=28)
plt.figtext(0.621, 0.36, 'G',clip_on=False,color='black', weight='bold',size=28)
###############################################################
gssub0 = gridspec.GridSpecFromSubplotSpec(1, 2 , subplot_spec=gs[0],wspace=0.2,hspace=0.05,width_ratios=[1,2]) #,height_ratios=[0.2,1])
gssub01 = gridspec.GridSpecFromSubplotSpec(2, 1 , subplot_spec=gssub0[0],wspace=0.05,hspace=0.1)

#ax1 = plt.subplot(gssub01[0])
#ax1.imshow(imBD)
#layoutOfPanel(ax1,xyInvisible=[True,True])

ax1 = plt.subplot(gssub01[1])

imMaskBD = np.zeros((opsBD['Ly'], opsBD['Lx']))

ax1.imshow(imBD,cmap='gray')

#intersectionROIs.append([n,m,xpixBD,ypixBD,xpixADPrime2,ypixADPrime2,intersection,eitherOr,intersection/eitherOr])
for n in ROIIDs:
    print(analysisBAndAD[n][0][0],analysisBAndAD[n][0][1])
    for i in intersectionROIs:
        #pdb.set_trace()
        if (i[0] == analysisBAndAD[n][0][0]) and (i[1] == analysisBAndAD[n][0][1]):
            print('match :',i[0],i[1])
            ypixBD = i[3]- aS.cutOffXPub
            xpixBD = i[2]- aS.cutOffXPub
            ax1.plot(xpixBD,ypixBD,zorder=1,alpha=0.6)



#ax1.imshow(imMaskBD,cmap='gist_ncar')

ax1.plot([265,390],[400,400],lw=8,c='w')
layoutOfPanel(ax1,xyInvisible=[True,True])
text11 = ax1.annotate(r'100 $\mu$m', xy=(295,435), annotation_clip=False,
           xytext=None, textcoords='data',fontsize=18,color='k',
           arrowprops=None
           )
###############################################################
# plotting example ROI traces and wheel activity
gssub02 = gridspec.GridSpecFromSubplotSpec(2, 1 , subplot_spec=gssub0[1],wspace=0.05,hspace=0.2)
gssub021 = gridspec.GridSpecFromSubplotSpec(2, len(recordingsBD) , subplot_spec=gssub02[0],wspace=0.05,hspace=0.1,height_ratios=[0.2,1])

figsBD = []
for n in range(len(recordingsBD)):
    ax0 = plt.subplot(gssub021[n])
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.spines['bottom'].set_visible(False)
    ax0.xaxis.set_visible(False)
    majorLocator_x = plt.MultipleLocator(100)
    ax0.yaxis.set_major_locator(majorLocator_x)
    if n == 0:
        ax0.spines['left'].set_position(('outward', 10))
        ax0.yaxis.set_ticks_position('left')
        ax0.set_ylabel('cm')
        #ax0.set_xlabel('time (s)')
    else:
        ax0.spines['left'].set_visible(False)
        ax0.yaxis.set_visible(False)  # set_ticks_position('left')


    ax1 = plt.subplot(gssub021[n+len(recordingsBD)])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_position(('outward', 10))
    ax1.xaxis.set_ticks_position('bottom')
    #ax1.set_xlabel('time (s)')
    if n == 0:
        ax1.spines['left'].set_position(('outward', 10))
        ax1.yaxis.set_ticks_position('left')
        ax1.set_ylabel(r'$\Delta F/F_{\rm rest}$')
    else:
        ax1.spines['left'].set_visible(False)
        ax1.yaxis.set_visible(False)  # set_ticks_position('left')
    ax0.set_ylim(-4,260)
    ax1.set_ylim(-3,40)
    text11 = ax0.annotate('rec. %s' % (recordingsBD[n]+1), xy=(0, 220), annotation_clip=False, xytext=None, textcoords='data', fontsize=14, color='k', arrowprops=None)
    figsBD.append([ax0,ax1])

nROI = 0
for i in ROIIDs: #(10): #len(analysisBAndAD)):
    nFig = 0
    for n in recordingsBD:
        #pdb.set_trace()
        if i==0:
            figsBD[nFig][0].axvspan(aS.activityTimes[0], aS.activityTimes[1], alpha=0.2, color='purple')
            figsBD[nFig][0].plot(wheelActivityBD[n][1][4],wheelActivityBD[n][1][3]* 80. / 360.,lw=2,c='k')
            figsBD[nFig][1].axvspan(aS.activityTimes[0], aS.activityTimes[1], alpha=0.2, color='purple')
        smoothedActivity = movingAverage(analysisBAndAD[i][1][n][4][:,2], 7)
        #figsBD[nFig][1].plot(analysisBAndAD[i][1][n][4][:,0],analysisBAndAD[i][1][n][4][:,2]+nROI*5,lw=1)
        figsBD[nFig][1].plot(analysisBAndAD[i][1][n][4][:,0],smoothedActivity + nROI*6, lw=2,alpha=0.7,clip_on=False)
        nFig+=1
    nROI+=1

gssub022 = gridspec.GridSpecFromSubplotSpec(2, len(recordingsAD) , subplot_spec=gssub02[1],wspace=0.05,hspace=0.1,height_ratios=[0.2,1]) #####################


# creating of figure instances
figsAD = []
for n in range(len(recordingsAD)):
    ax0 = plt.subplot(gssub022[n])
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.spines['bottom'].set_visible(False)
    ax0.xaxis.set_visible(False)
    majorLocator_x = plt.MultipleLocator(100)
    ax0.yaxis.set_major_locator(majorLocator_x)
    if n == 0:
        ax0.spines['left'].set_position(('outward', 10))
        ax0.yaxis.set_ticks_position('left')
        ax0.set_ylabel('cm')
    else:
        ax0.spines['left'].set_visible(False)
        ax0.yaxis.set_visible(False)  # set_ticks_position('left')

    ax1 = plt.subplot(gssub022[n+len(recordingsAD)])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_position(('outward', 10))
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_xlabel('Time (s)')
    if n == 0:
        ax1.spines['left'].set_position(('outward', 10))
        ax1.yaxis.set_ticks_position('left')
        ax1.set_ylabel(r'$\Delta F/F_{\rm rest}$')
        #ax1.set_xlabel('time (s)')
    else:
        ax1.spines['left'].set_visible(False)
        ax1.yaxis.set_visible(False)  # set_ticks_position('left')
    ax0.set_ylim(-4,260)
    ax1.set_ylim(-3, 40)
    text11 = ax0.annotate('rec. %s' % (trialsBeforeDrug+recordingsAD[n]+1), xy=(0, 220), annotation_clip=False, xytext=None, textcoords='data', fontsize=14, color='k', arrowprops=None)
    figsAD.append([ax0,ax1])

nROI = 0
for i in ROIIDs: #(10): #len(analysisBAndAD)):
    nFig = 0
    for n in recordingsAD: #(len(analysisBAndAD[i][2])):
        #pdb.set_trace()
        if nROI==0:
            figsAD[nFig][0].plot(wheelActivityAD[n][1][4],wheelActivityAD[n][1][3]* 80. / 360.,lw=2,c='k')
            figsAD[nFig][0].axvspan(aS.activityTimes[0], aS.activityTimes[1], alpha=0.2, color='purple')
            figsAD[nFig][1].axvspan(aS.activityTimes[0], aS.activityTimes[1], alpha=0.2, color='purple')
        smoothedActivity = movingAverage(analysisBAndAD[i][2][n][4][:,2], 7)
        #figsAD[nFig][1].plot(analysisBAndAD[i][2][n][4][:,0],analysisBAndAD[i][2][n][4][:,2]+nROI*5,lw=1)
        figsAD[nFig][1].plot(analysisBAndAD[i][2][n][4][:, 0],smoothedActivity + nROI * 6, lw=2,clip_on=False)
        nFig+=1
    nROI+=1

###############################################################
# plotting mean activity vs. trial number
gssub1 = gridspec.GridSpecFromSubplotSpec(1, 3 , subplot_spec=gs[1],wspace=0.2,hspace=0.1,width_ratios=[1,1,1.2]) #,height_ratios=[0.2,1])

ax0 = plt.subplot(gssub1[0])
#ax0.set_title('mean activity during locomotion')
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_position(('outward', 10))
ax0.xaxis.set_ticks_position('bottom')
ax0.spines['left'].set_position(('outward', 10))
ax0.yaxis.set_ticks_position('left')
ax0.set_xlabel('Recording number')
ax0.set_ylabel(r'Locomotion fluorescence ($\Delta F/F_{\rm rest}$)')
majorLocator_x = plt.MultipleLocator(1)
ax0.xaxis.set_major_locator(majorLocator_x)
bracket = ax0.annotate('', xy=(1, 9),  xycoords='data', annotation_clip=False,
           xytext=(6, 9), textcoords='data',
           arrowprops=dict(arrowstyle="-",
                           connectionstyle='bar,fraction=0.2', linewidth=2,
                           ec='k',
                           shrinkA=10, shrinkB=10,
                           )
           )
bracket = ax0.annotate('', xy=(7, 9),  xycoords='data', annotation_clip=False,
           xytext=(14, 9), textcoords='data',
           arrowprops=dict(arrowstyle="-",
                           connectionstyle='bar,fraction=0.14', linewidth=2,
                           ec='k',
                           shrinkA=10, shrinkB=10,
                           )
           )
text11 = ax0.annotate('before drug', xy=(1.3,9.8), annotation_clip=False,
           xytext=None, textcoords='data',fontsize=18,color='k',
           arrowprops=None
           )
text11 = ax0.annotate('after drug', xy=(8.5,9.8), annotation_clip=False,
           xytext=None, textcoords='data',fontsize=18,color='k',
           arrowprops=None
           )

ax1 = plt.subplot(gssub1[1])
#ax1.set_title('mean activity during locomotion')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.xaxis.set_ticks_position('bottom')
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.set_xlabel('Recording number')
ax1.set_ylabel(r'Normalized locomotion fluorescence ($\Delta F/F_{\rm rest}$)')
ax1.xaxis.set_major_locator(majorLocator_x)


mColors = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9']
fit = 'exponential'
ni=0
for i in ROIIDs:
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

    if ni==0:
        ax0.axhline(y=0,ls='--',c='0.5')
        #ax1.axhline(y=0, ls='--', c='0.5')

    #if i<10:
    ax0.plot(range(1,trialsBeforeDrug+trialsAfterDrug+1),activity,'o-',ms=4,c=mColors[ni],alpha=0.3,lw=2)
    ax0.plot(range(trialsBeforeDrug + 1, trialsBeforeDrug + trialsAfterDrug + 1), activity[trialsBeforeDrug:], 'o-', ms=4,c=mColors[ni],lw=2)

    if ni==0:
        ax1.axhline(y=0,ls='--',c='0.7')
        #ax1.axhline(y=0, ls='--', c='0.7')

    ax1.plot(range(1, trialsBeforeDrug + trialsAfterDrug + 1), activity - activityFit, 'o-', ms=4, c=mColors[ni], alpha=0.3,lw=2)

    ax1.plot(range(trialsBeforeDrug + 1, trialsBeforeDrug + trialsAfterDrug + 1), activity[trialsBeforeDrug:] - activityFit[trialsBeforeDrug:], 'o-', ms=4, c=mColors[ni],lw=2)

    ni+=1

########################################################################################
# summary of all animals

#allData = [trialsBeforeDrug,trialsAfterDrug,integral,activity,integralEv,activityEv,F0Ev]
#animals = ['animal#1','animal#3','animal#2','animal#4','animal#1_2']
animals = ['animal#3','animal#2','animal#4','animal#1_2']

allData = []
for a in range(len(animals)):
    aD = pickle.load(open( dataOutDir + 'activityChanges_%s.p' % animals[a], 'rb' ) )
    allData.append(aD)



ax1 = plt.subplot(gssub1[2]) ############################################################
#a = [0,3]
#ax1.plot(a,a,c='0.5')
symbols = ['o','v','>','s','o']
ax1.axhline(y=0,ls='--',c='0.5')
cum = []
nAnimal = 0
nROIs = 0
#for n in range(len(allData)):
for a in range(len(animals)):
    allData = pickle.load(open( dataOutDir + 'activityChanges_%s.p' % animals[a], 'rb' ) )
    for i in range(len(allData[5])):
        #pdb.set_trace()
        if allData[5][i,2]<0.01 and allData[5][i,3] == 1:
            ax1.plot(allData[5][i,0],allData[5][i,1],symbols[nAnimal],ms=7,c='blue')
        elif allData[5][i,2]<0.01 and allData[5][i,3] == 0:
            ax1.plot(allData[5][i,0],allData[5][i,1],symbols[nAnimal],ms=7,c='red')
        else:
            ax1.plot(allData[5][i,0],allData[5][i,1],symbols[nAnimal],c='gray',ms=7,alpha=0.5)
        cum.append((allData[5][i,0],allData[5][i,1]))
        nROIs+=1
    nAnimal +=1

print('total number of ROIs : ', nROIs)
cum = np.asarray(cum)
mask = cum[:,0]<10.
out = linregress(cum[:,0][mask], cum[:,1][mask])
meanBefore = cum[:,0][mask]
ax1.plot(meanBefore,out[1]+out[0]*meanBefore,c='0.6')

print(out)
text11 = ax1.annotate(r'$r= %s$, R$^2$ = %s, p < 0.0001' %(np.round(out[2],3),np.round(out[2]**2,3)) , xy=(0,-6), annotation_clip=False,
           xytext=None, textcoords='data',fontsize=16,
           arrowprops=None
           )
#ax1.plot(activityEv[:,0],activityEv[:,1],'o',ms=2,alpha=0.5)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.xaxis.set_ticks_position('bottom')
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.set_xlabel(r'Locomotion fluorescence before drug ($\Delta F/F_{\rm rest}$)')
ax1.set_ylabel(r'Change in locomotion fluo. after drug ($\Delta F/F_{\rm rest}$)')
#ax1.set_xlim(-2,10)
#ax1.set_ylim(-10,5)


########################################################################################

plt.savefig(figOutDir+'WalkingFluorescenceBeforeAfterDrug_2.pdf')
plt.savefig(figOutDir+'WalkingFluorescenceBeforeAfterDrug_2.svg')
