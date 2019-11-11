import numpy as np
import matplotlib.pyplot as plt

#########################################################
# set parameters

#animalID = 'animal#1_2' #'animal#3' #'animal#1' #'animal#3' # 'animal#2'


class animalSettings:
    def __init__(self, mouse):
        animals = ['animal#1','animal#3','animal#2','animal#4','animal#1_2']
        if mouse in animals:
            self.animalID = mouse
        else:
            print('Mouse %s does not exist.' % mouse)

    def loadSettings(self):

        if self.animalID == 'animal#1':

            self.baseDir = '/media/labo_rw/JINmGluR/'
            self.beforeDrugDir = '2019.08.15_002_animal#1/animal#1_00000-00013/'
            self.afterDrugDir = '2019.08.15_002_animal#1/animal#1_00015-20_23-28/'
            self.suite2pDir = 'suite2p/plane0/'

            self.cutOffX = 8
            self.cutOffY = 19

            self.xOffset =-30.
            self.yOffset = 6.

            self.baselinePeriodsBD = [[0, 13]]
            self.baselinePeriodsAD = [[0, 11]]

            self.nImgsBD = [[0, '000-000'], [1, '000-001'], [2, '000-002'], [3, '000-003'], [4, '000-004'], [5, '000-005'], [6, '001-000'], [7, '001-001'], [8, '002-000'], [9, '002-001'],
                       [10, '002-002'], [11, '002-003'], [12, '002-004'], [13, '002-005']]  # range(startBD,endBD+1)
            self.nImgsAD = [[15, '004-000'], [16, '004-001'], [17, '004-002'], [18, '004-003'], [19, '004-004'], [20, '004-005'], [23, '007-000'], [24, '007-001'], [25, '007-002'], [26, '007-003'],
                       [27, '007-004'], [28, '007-005']]

            self.minOverlap = 0.5
            self.intLims = None

        elif self.animalID == 'animal#3':
            self.baseDir = '/media/labo_rw/JINmGluR/'
            self.beforeDrugDir = '2019.08.16_000_animal#3/animal#3_00001-6_8-13/'
            self.afterDrugDir = '2019.08.16_000_animal#3/animal#3_00018_20_22_24_26_28_30_32_34/'
            self.suite2pDir = 'suite2p/plane0/'

            self.cutOffX = 49
            self.cutOffY = 32

            self.xOffset = 0
            self.yOffset = 0

            self.baselinePeriodsBD = [[0,5],[6,11]]
            self.baselinePeriodsAD = [[0,8]]
            self.nImgsBD = [[1, '001-000'], [2, '001-001'], [3, '001-002'], [4, '001-003'], [5, '001-004'], [6, '001-005'],[8, '003-000'], [9, '003-001'], [10, '003-002'], [11, '003-003'], [12, '003-004'], [13, '003-005']]  # range(startBD,endBD+1)
            self.nImgsAD = [[18, 'Jin-004'], [20, 'Jin-006'], [22, 'Jin-007'], [24, 'Jin-008'], [26, 'Jin-009'], [28, 'Jin-010'], [30, 'Jin-011'], [32, 'Jin-012'], [34, 'Jin-013']]

            self.minOverlap = 0.4 # fraction of minial ROI overlap to consider the ROIs as identical
            self.intLims = None

        elif self.animalID == 'animal#2':
            self.baseDir = '/media/labo_rw/JINmGluR/'
            self.beforeDrugDir = '2019.08.16_001_animal#2/animal#2_00000-5/'
            self.afterDrugDir = '2019.08.16_001_animal#2/animal#2_00014_16_18_20_22_24_26_28/'
            self.suite2pDir = 'suite2p/plane0/'

            self.cutOffX = 49
            self.cutOffY = 32

            self.xOffset = 0
            self.yOffset = 0

            self.baselinePeriodsBD = [[0,5]]
            self.baselinePeriodsAD = [[0,7]]
            self.nImgsBD = [[0, '000-000'], [1, '000-001'], [2, '000-002'], [3, '000-003'], [4, '000-004'], [5, '000-005']] # range(startBD,endBD+1)
            self.nImgsAD = [[14, 'Jin-002'], [16, 'Jin-003'], [18, 'Jin-004'], [20, 'Jin-005'], [22, 'Jin-006'], [24, 'Jin-007'], [26, 'Jin-008'], [28, 'Jin-009']]

            self.minOverlap = 0.3 # fraction of minial ROI overlap to consider the ROIs as identical

            self.intLims = [-60,40,-10,150]


        elif self.animalID == 'animal#4':
            self.baseDir = '/media/labo_rw/JINmGluR/'
            self.beforeDrugDir = '2019.08.17_000_animal#4/animal#4_00000-5/'
            self.afterDrugDir = '2019.08.17_000_animal#4/animal#4_00014_16_18_20_22_24_26_28/'
            self.suite2pDir = 'suite2p/plane0/'

            self.cutOffX = 49
            self.cutOffY = 50

            self.xOffset = 0
            self.yOffset = 0

            self.baselinePeriodsBD = [[0,5]]
            self.baselinePeriodsAD = [[0,7]]
            self.nImgsBD = [[0, '000-000'], [1, '000-001'], [2, '000-002'], [3, '000-003'], [4, '000-004'], [5, '000-005']] # range(startBD,endBD+1)
            self.nImgsAD = [[14, 'Jin-002'], [16, 'Jin-003'], [18, 'Jin-004'], [20, 'Jin-005'], [22, 'Jin-006'], [24, 'Jin-007'], [26, 'Jin-008'], [28, 'Jin-009']]

            self.minOverlap = 0.3 # fraction of minial ROI overlap to consider the ROIs as identical
            self.intLims = None

        elif self.animalID == 'animal#1_2':
            self.baseDir = '/media/labo_rw/JINmGluR/'
            self.beforeDrugDir = '2019.08.17_002_animal#1/animal#1_00000-5/'
            self.afterDrugDir = '2019.08.17_002_animal#1/animal#1_00014_16_18_20_22_24_26/'
            self.suite2pDir = 'suite2p/plane0/'

            self.cutOffX = 60
            self.cutOffY = 53

            self.xOffset = -38
            self.yOffset = 10

            self.baselinePeriodsBD = [[0,5]]
            self.baselinePeriodsAD = [[0,7]]
            self.nImgsBD = [[0, '000-000'], [1, '000-001'], [2, '000-002'], [3, '000-003'], [4, '000-004'], [5, '000-005']] # range(startBD,endBD+1)
            self.nImgsAD = [[14, 'Jin-002'], [16, 'Jin-003'], [18, 'Jin-004'], [20, 'Jin-005'], [22, 'Jin-006'], [24, 'Jin-007'], [26, 'Jin-008']]

            self.minOverlap = 0.4 # fraction of minial ROI overlap to consider the ROIs as identical
            self.intLims = None
