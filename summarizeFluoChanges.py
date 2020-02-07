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
#from sklearn.linear_model import LinearRegression
from scipy.stats import linregress
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

#allData = [trialsBeforeDrug,trialsAfterDrug,integral,activity,integralEv,activityEv,F0Ev]
#animals = ['animal#1','animal#3','animal#2','animal#4','animal#1_2']
animals = ['animal#3','animal#2','animal#4','animal#1_2']

allData = []
for a in range(len(animals)):
    aD = pickle.load(open( dataOutDir + 'activityChanges_%s.p' % animals[a], 'rb' ) )
    allData.append(aD)


##################################################################
#GT.generateSingleTraceFigure(aS.animalID,timeStampsBD,timeStampsAD,analysisBAndAD,aS.baselinePeriodsBD,aS.baselinePeriodsAD,columnIdx=1)
#GT.generateSingleTraceFigure(aS.animalID,timeStampsBD,timeStampsAD,analysisBAndAD,aS.baselinePeriodsBD,aS.baselinePeriodsAD,columnIdx=2)
#GT.generateSingleTraceFigure(aS.animalID,timeStampsBD,timeStampsAD,analysisBAndAD,aS.baselinePeriodsBD,aS.baselinePeriodsAD,columnIdx=3)
#pdb.set_trace()
##################################################################
# Show final results

fig_width = 10 # width in inches
fig_height = 4  # height in inches
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
gs = gridspec.GridSpec(1, 2,
                       width_ratios=[1.5,1],
                       #height_ratios=[1,1,2.5,1]
                       )
# define vertical and horizontal spacing between panels
gs.update(wspace=0.3,hspace=0.25)

# possibly change outer margins of the figure
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

plt.figtext(0.06, 0.95, 'summary of fluorescence changes across all animals', clip_on=False, color='black', size=12)


ax1 = plt.subplot(gs[0]) ############################################################
#a = [0,3]
#ax1.plot(a,a,c='0.5')
symbols = ['o','v','>','s','D']
ax1.axhline(y=0,ls='--',c='0.5')
cum = []
nAnimal = 0
#for n in range(len(allData)):
for a in range(len(animals)):
    allData = pickle.load(open( dataOutDir + 'activityChanges_%s.p' % animals[a], 'rb' ) )
    for i in range(len(allData[5])):
        #pdb.set_trace()
        if allData[5][i,2]<0.01 and allData[5][i,3] == 1:
            ax1.plot(allData[5][i,0],allData[5][i,1],symbols[nAnimal],ms=3,c='C0')
        elif allData[5][i,2]<0.01 and allData[5][i,3] == 0:
            ax1.plot(allData[5][i,0],allData[5][i,1],symbols[nAnimal],ms=3,c='C1')
        else:
            ax1.plot(allData[5][i,0],allData[5][i,1],symbols[nAnimal],c='C2',ms=2,alpha=0.5)
        cum.append((allData[5][i,0],allData[5][i,1]))
    nAnimal +=1


cum = np.asarray(cum)
mask = cum[:,0]<10.
out = linregress(cum[:,0][mask], cum[:,1][mask])
meanBefore = cum[:,0][mask]
ax1.plot(meanBefore,out[1]+out[0]*meanBefore,c='0.6')
print(out)
text11 = ax1.annotate(r'R$^2$ = %s, p < 0.0001' %(np.round(out[2]**2,3)) , xy=(-1,-9), annotation_clip=False,
           xytext=None, textcoords='data',fontsize=10,
           arrowprops=None
           )
#ax1.plot(activityEv[:,0],activityEv[:,1],'o',ms=2,alpha=0.5)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.xaxis.set_ticks_position('bottom')
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.set_xlabel('mean activity before drug during locomotion')
ax1.set_ylabel('change in activity after drug during locomotion')
#ax1.set_xlim(-2,10)
#ax1.set_ylim(-10,5)



ax1 = plt.subplot(gs[1]) ############################################################
#a = [0,3]
#ax1.plot(a,a,c='0.5')
ax1.axhline(y=0,ls='--',c='0.5')
allPairs = []
nAnimal = 0
#for n in range(len(allData)):
for a in range(len(animals)):
    allData = pickle.load(open( dataOutDir + 'activityChanges_%s.p' % animals[a], 'rb' ) )
    for i in range(len(allData[5])):
        #pdb.set_trace()
        ax1.plot([1,2], [allData[5][i, 4],allData[5][i, 5]], '%s-' % symbols[nAnimal],lw=0.5, ms=3, c='C0',alpha=0.2)
        allPairs.append([allData[5][i, 4],allData[5][i, 5]])
    nAnimal+=1

allPairs = np.asarray(allPairs)
ptt  = stats.ttest_rel(allPairs[:,0],allPairs[:,1])
print('paired t-test :',ptt)
mask = allPairs[:,0] > 1
ptt2  = stats.ttest_rel(allPairs[:,0][mask],allPairs[:,1][mask])
print('paired t-test large values (>1) :',ptt2)

text11 = ax1.annotate('paired t-test: p = %s\npaired t-test (activity before>1): p = %s' % (np.round(ptt[1],4),np.round(ptt2[1],4)) , xy=(1,-1.1), annotation_clip=False,
           xytext=None, textcoords='data',fontsize=10,
           arrowprops=None
           )
#ax1.plot(activityEv[:,0],activityEv[:,1],'o',ms=2,alpha=0.5)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.xaxis.set_ticks_position('bottom')
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.set_ylabel('mean activity')
majorLocator_x = plt.MultipleLocator(1)
ax1.xaxis.set_major_locator(majorLocator_x)
labels = [item.get_text() for item in ax1.get_xticklabels()]
labels[1] = 'before'
labels[2] = 'after'
ax1.set_xticklabels(labels)

#ax1.set_ylabel('before vs after')
#ax1.set_xlim(-2,10)
#ax1.set_ylim(-10,5)



plt.savefig(figOutDir+'SummaryFluorescenceChanges.pdf')
#plt.savefig(figOutDir+'FluorescenceTraces_%s.png' % animalID)
#plt.show()



