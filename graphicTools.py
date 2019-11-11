import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
import pdb

#########################################################
# set parameters

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


def generateSingleTraceFigure(animalID,timeStampsBD,timeStampsAD,analysisBAndAD,periodBD,periodAD,columnIdx = 1):
    figOutDir = 'figureOutput/'

    nRecsBD = len(timeStampsBD)
    nRecsAD = len(timeStampsAD)
    #print('nRecsAD',nRecsAD)

    periods = []
    for i in range(nRecsBD):
        pN = 0
        for per in periodBD:
            if per[0]<=i and i<=per[1]:
                periods.append([i,pN])
            pN+=1
    for i in range(nRecsAD):
        pNAD = pN
        for per in periodAD:
            if per[0]<=i and i<=per[1]:
                periods.append([i+nRecsBD,pNAD])
            pNAD+=1

    #pdb.set_trace()

    #####################################################################
    fig_width = 50  # width in inches
    fig_height = len(analysisBAndAD)*3  # height in inches
    fig_size = [fig_width, fig_height]
    params = {'axes.labelsize': 11, 'axes.titlesize': 11, 'font.size': 11, 'xtick.labelsize': 11, 'ytick.labelsize': 11, 'figure.figsize': fig_size, 'savefig.dpi': 600,
              'axes.linewidth': 1.3, 'ytick.major.size': 4,  # major tick size in points
              'xtick.major.size': 4  # major tick size in points
              # 'edgecolor' : None
              # 'xtick.major.size' : 2,
              # 'ytick.major.size' : 2,
              }
    rcParams.update(params)
    # create figure instance
    fig = plt.figure()

    # define sub-panel grid and possibly width and height ratios
    gs = gridspec.GridSpec(len(analysisBAndAD), nRecsBD+nRecsAD,  # width_ratios=[1,1.2],
                           #height_ratios=[1, 1, 2.5, 1]
                           )
    # define vertical and horizontal spacing between panels
    gs.update(wspace=0.3, hspace=0.25)

    # possibly change outer margins of the figure
    plt.subplots_adjust(left=0.05, right=0.98, top=0.98, bottom=0.02)
    plt.figtext(0.06, 0.99, '%s : %s ROIs,   %s recs. w/o drug; %s recordings with drug' % (animalID, len(analysisBAndAD), len(analysisBAndAD[0][1]), len(analysisBAndAD[0][2])), clip_on=False, color='black', size=14)

    #gssub0 = gridspec.GridSpecFromSubplotSpec(len(analysisBAndAD), nRecsBD+nRecsAD, subplot_spec=gs[0], wspace=0.05, hspace=0.1, height_ratios=[0.2, 1])

    # creating of figure instances
    allFigs = []
    roiRange = []
    nFigs = 0

    for i in range(len(analysisBAndAD)):
        figsPerRoi = []
        roiMinMax = [100,-100]
        for n in range(nRecsBD):
            ax0 = plt.subplot(gs[n+nFigs])
            layoutOfPanel(ax0)
            figsPerRoi.append(ax0)
            if min(analysisBAndAD[i][1][n][4][:, columnIdx])<roiMinMax[0]:
                roiMinMax[0] = min(analysisBAndAD[i][1][n][4][:, columnIdx])
            if max(analysisBAndAD[i][1][n][4][:, columnIdx])>roiMinMax[1]:
                roiMinMax[1] = max(analysisBAndAD[i][1][n][4][:, columnIdx])
        for n in range(nRecsAD):
            ax1 = plt.subplot(gs[n+nRecsBD+nFigs])
            layoutOfPanel(ax1)
            figsPerRoi.append(ax1)
            #print(i,n)
            if min(analysisBAndAD[i][2][n][4][:, columnIdx])<roiMinMax[0]:
                roiMinMax[0] = min(analysisBAndAD[i][2][n][4][:, columnIdx])
            if max(analysisBAndAD[i][2][n][4][:, columnIdx])>roiMinMax[1]:
                roiMinMax[1] = max(analysisBAndAD[i][2][n][4][:, columnIdx])
        roiRange.append(roiMinMax)

        nFigs += (nRecsBD+nRecsAD)
        allFigs.append(figsPerRoi)
    ccc = ['C0','C1','C2','C3']
    for i in range(len(analysisBAndAD)):
        for n in range(nRecsBD):
            allFigs[i][n].axhline(y=0,ls='--',c='0.5')
            allFigs[i][n].plot(analysisBAndAD[i][1][n][4][:, 0], analysisBAndAD[i][1][n][4][:, columnIdx],c=ccc[periods[n][1]])
            allFigs[i][n].set_ylim(roiRange[i][0],roiRange[i][1])
        for n in range(nRecsAD):
            allFigs[i][nRecsBD+n].axhline(y=0, ls='--', c='0.5')
            allFigs[i][nRecsBD+n].plot(analysisBAndAD[i][2][n][4][:, 0], analysisBAndAD[i][2][n][4][:, columnIdx],c=ccc[periods[n+nRecsBD][1]])
            allFigs[i][nRecsBD+n].set_ylim(roiRange[i][0], roiRange[i][1])

            #figsAD[n][1].plot(analysisBAndAD[i][2][n][4][:,0],analysisBAndAD[i][2][n][4][:,2]+i*2,lw=0.2)
    columnIdentity = {1:'rawFluo',2:'deltaFOverF0',3:'deltaF'}
    plt.savefig(figOutDir + 'AllFluorescenceTraces_%s_%s.pdf' % (animalID,columnIdentity[columnIdx]))
    # plt.savefig(figOutDir+'FluorescenceTraces_%s.png' % animalID)
    #plt.show()