import numpy as np
import matplotlib.pyplot as plt
import cv2


stat = np.load('stat.npy')
ops = np.load('ops.npy').item()

im = np.zeros((ops['Ly'], ops['Lx']))

for n in range(0,ncells):
    ypix = stat[n]['ypix'][~stat[n]['overlap']]
    xpix = stat[n]['xpix'][~stat[n]['overlap']]
    im[ypix,xpix] = n+1