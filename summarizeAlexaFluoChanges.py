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
pubOutDir = 'publicationFigures/'
# animalID = 'animal#1'
# #baseDir = '/media/labo_rw/JINmGluR/'
# baseDir = '/home/mgraupe/'
# beforeDrugDir = '2019.08.15_002_animal#1/animal#1_00000-00013/'
# afterDrugDir = '2019.08.15_002_animal#1/animal#1_00015-20_23-28/'
# suite2pDir = 'suite2p/plane0/'

##########################################################
# read data and determine principal parameters

#allData = [trialsBeforeDrug,trialsAfterDrug,integral,activity,integralEv,activityEv,F0Ev]
animals = ['animal#1','animal#1_2','animal#2','animal#3','animal#4']
#animals = ['animal#3','animal#2','animal#4','animal#1_2']

allData = []
fitData = np.zeros(2)
for a in range(len(animals)):
    aD = pickle.load(open( dataOutDir + 'analysisOfAlexaRecordings_%s.p' % animals[a], 'rb' ) )
    allData.append([a,animals[a],aD])
    if animals[a] != 'animal#1':
        fitData = np.vstack((fitData,np.column_stack((aD[:,1],aD[:,2]))))
fitData = fitData[1:]

mask = fitData[:,0]>0.
maskN = fitData[:,0]<0.
B = np.mean(fitData[:,1][maskN])

fitfunc = lambda p, x: p[0]*(1. - np.exp(-x/p[1])) + B
errfunc = lambda p, x, y: fitfunc(p, x) - y

p0Int = [6.,1.]
p1Int, success = scipy.optimize.leastsq(errfunc, p0Int, args=(fitData[:,0][mask], fitData[:,1][mask]))
print(p1Int)
tt = np.linspace(0,40,100)
yy = fitfunc(p1Int,tt)

##################################################################
#GT.generateSingleTraceFigure(aS.animalID,timeStampsBD,timeStampsAD,analysisBAndAD,aS.baselinePeriodsBD,aS.baselinePeriodsAD,columnIdx=1)
#GT.generateSingleTraceFigure(aS.animalID,timeStampsBD,timeStampsAD,analysisBAndAD,aS.baselinePeriodsBD,aS.baselinePeriodsAD,columnIdx=2)
#GT.generateSingleTraceFigure(aS.animalID,timeStampsBD,timeStampsAD,analysisBAndAD,aS.baselinePeriodsBD,aS.baselinePeriodsAD,columnIdx=3)
#pdb.set_trace()
##################################################################
# Show final results

fig_width = 7 # width in inches
fig_height = 5  # height in inches
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
gs = gridspec.GridSpec(1,1,
                       #width_ratios=[1.5,1],
                       #height_ratios=[1,1,2.5,1]
                       )
# define vertical and horizontal spacing between panels
gs.update(wspace=0.3,hspace=0.25)

# possibly change outer margins of the figure
plt.subplots_adjust(left=0.11, right=0.95, top=0.95, bottom=0.15)

#plt.figtext(0.06, 0.95, 'summary of Alexa 594 fluorescence changes across all animals', clip_on=False, color='black', size=12)

animalLabels = ['animal#1A','animal#1','animal#2','animal#3','animal#4']
ax1 = plt.subplot(gs[0]) ############################################################
#a = [0,3]
#ax1.plot(a,a,c='0.5')
symbols = ['o','v','>','s','D']
ax1.axvline(x=0,ls='--',c='0.5')
cc = ['0.6','C0','C1','C2','C3']
for n in range(len(allData)):
    ax1.errorbar(allData[n][2][:,1],allData[n][2][:,2],allData[n][2][:,3],fmt='o-',c=cc[n],ms=2,label='%s' % animalLabels[n])

ax1.plot(tt,yy,c='k',lw=3,label=r'exp. fit : $\tau = %s$ min' % np.round(p1Int[1],1))
#ax1.plot(activityEv[:,0],activityEv[:,1],'o',ms=2,alpha=0.5)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.xaxis.set_ticks_position('bottom')
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')

plt.legend(loc=(0.02,0.65))
# change legend text size
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=10)

ax1.set_xlabel('Time since Alex 594 application (min)')
ax1.set_ylabel('Raw fluorescence')
#ax1.set_xlim(-2,10)
#ax1.set_ylim(-10,5)




plt.savefig(pubOutDir+'SummaryAlexa594FluorescenceChanges_publication.pdf')
#plt.savefig(figOutDir+'FluorescenceTraces_%s.png' % animalID)
#plt.show()



