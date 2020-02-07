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
animals = ['animal#1','animal#3','animal#2','animal#4','animal#1_2']
#animals = ['animal#3','animal#2','animal#4','animal#1_2']

allData = []
for a in range(len(animals)):
    aD = pickle.load(open( dataOutDir + 'analysisBeforeAndAt820_%s.p' % animals[a], 'rb' ) )
    allData.append(aD)


##################################################################
#GT.generateSingleTraceFigure(aS.animalID,timeStampsBD,timeStampsAD,analysisBAndAD,aS.baselinePeriodsBD,aS.baselinePeriodsAD,columnIdx=1)
#GT.generateSingleTraceFigure(aS.animalID,timeStampsBD,timeStampsAD,analysisBAndAD,aS.baselinePeriodsBD,aS.baselinePeriodsAD,columnIdx=2)
#GT.generateSingleTraceFigure(aS.animalID,timeStampsBD,timeStampsAD,analysisBAndAD,aS.baselinePeriodsBD,aS.baselinePeriodsAD,columnIdx=3)
#pdb.set_trace()
##################################################################
# Show final results

fig_width = 6 # width in inches
fig_height = 5  # height in inches
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
gs = gridspec.GridSpec(1, 1,
                       #width_ratios=[1.5,1],
                       #height_ratios=[1,1,2.5,1]
                       )
# define vertical and horizontal spacing between panels
gs.update(wspace=0.3,hspace=0.25)

# possibly change outer margins of the figure
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

plt.figtext(0.06, 0.95, 'summary of fluctuations during 910 and 820 nm runs', clip_on=False, color='black', size=12)


ax1 = plt.subplot(gs[0]) ############################################################
#a = [0,3]
#ax1.plot(a,a,c='0.5')
symbols = ['o','v','>','s','D']
#ax1.axhline(y=0,ls='--',c='0.5')
#cum = []
#nAnimal = 0
#for n in range(len(allData)):
aa = []
for a in range(len(animals)): # loop over all animals
    ratioPerAnimal = []
    for b in range(len(allData[a])): # loop over all intersection ROIs
        stdBD = []
        std820= []
        print('%s pre-durg and %s recording at 820' %(len(allData[a][b][1]),len(allData[a][b][2])))
        for c in range(len(allData[a][b][1])): # loop over all pre-drug ROIs
            #pdb.set_trace()
            stdBD.append(allData[a][b][1][c][2])
        for d in range(len(allData[a][b][2])): # loop over all 820 ROIs
            #pdb.set_trace()
            std820.append(allData[a][b][2][d][2])
        ratio = np.mean(std820)/np.mean(stdBD)
        ratioPerAnimal.append(ratio)
    ratioPerAnimal = np.asarray(ratioPerAnimal)
    ratioPerAnimal = ratioPerAnimal[ratioPerAnimal<3]
    aa.append(ratioPerAnimal)
#print(a,b)
#pdb.set_trace()
ax1.boxplot(aa)
#ax1.plot(a+(np.random.rand(len(ratioPerAnimal))-0.5)*0.2,ratioPerAnimal,symbols[a],ms=3,c='0.5')
#ax1.plot(a,np.mean(ratioPerAnimal),symbols[a],ms=5,c='C0')

#ax1.plot(activityEv[:,0],activityEv[:,1],'o',ms=2,alpha=0.5)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.xaxis.set_ticks_position('bottom')
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.set_ylabel('STD ratio (820 nm/910 nm before drug)')
majorLocator_x = plt.MultipleLocator(1)
ax1.xaxis.set_major_locator(majorLocator_x)
labels = [item.get_text() for item in ax1.get_xticklabels()]
print(labels)
for i in range(len(animals)):
    print(i)
    labels[i+1] = animals[i]
#labels[1] = 'before'
#labels[2] = 'after'
ax1.set_xticklabels(labels)

#ax1.set_ylabel('before vs after')
#ax1.set_xlim(-2,10)
#ax1.set_ylim(-10,5)



plt.savefig(figOutDir+'Summary820Fluctuations.pdf')
#plt.savefig(figOutDir+'FluorescenceTraces_%s.png' % animalID)
#plt.show()



