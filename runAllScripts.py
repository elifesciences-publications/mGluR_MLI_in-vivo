import numpy as np
import matplotlib.pyplot as plt
import sys, os

animals = ['animal#1','animal#3','animal#2','animal#4','animal#1_2']


for a in animals:
    print('analyzing animal %s' % a)
    #os.system('python alignImages.py %s' % a)
    #os.system('python matchROIs.py %s' % a)
    os.system('python analyzeFluoTraces.py %s' % a)


