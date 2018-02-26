#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
import argparse
import sys
import itertools

parser = argparse.ArgumentParser(description='Plot Chi^2 Contours')
parser.add_argument('filename', metavar='File', type=str, help='File containing chi^2 data')
parser.add_argument('cols', metavar='cols', type=str, nargs='+', help = 'Name of Columns')

args = parser.parse_args()

combos = itertools.combinations(args.cols, 2)

dtList = [(c, np.float64) for c in args.cols]
dtList.append(('Chi^2', np.float64))
dt = np.dtype(dtList)

data = np.loadtxt(args.filename, dt)

minChisq = np.amin(data['Chi^2'])
maxChisq = minChisq + 4

minPoint = data[data['Chi^2'] == minChisq]

for (col1, col2) in combos:
    cutCols = list(args.cols) #got to copy!
    cutCols.remove(col1)
    cutCols.remove(col2)
    
    dataMask = data['Chi^2'] < minChisq + 4
    for c in cutCols:
        dataMask = np.logical_and(dataMask, data[c] == minPoint[0][c])
    
    minChisqBasin = data[dataMask]

    minCol1 = np.amin(minChisqBasin[col1])
    maxCol1 = np.amax(minChisqBasin[col1])
    minCol2 = np.amin(minChisqBasin[col2])
    maxCol2 = np.amax(minChisqBasin[col2])

    #interpRegion = data[data[col1] > minCol1 and data[col1] < maxCol1 and data[col2] > minCol2 and data[col2] < maxCol2]

    xs = np.linspace(minCol1, maxCol1, 200)
    ys = np.linspace(minCol2, maxCol2, 200)
    xGrid, yGrid = np.meshgrid(xs, ys)
    zsInterp = interp.griddata((minChisqBasin[col1], minChisqBasin[col2]), minChisqBasin['Chi^2'], (xGrid, yGrid), method='linear')
    
    contours = [minChisq + x for x in np.arange(0.0, 4.1, 0.2)]

    plt.clf()
    contour_set_bands = plt.contour(xGrid, yGrid, zsInterp-minChisq, [1, 2], colors='k')
    contour_set_full = plt.contourf(xGrid, yGrid, zsInterp, contours, norm=plt.Normalize(vmax=maxChisq, vmin=minChisq))
#    plt.scatter(minChisqBasin[col1], minChisqBasin[col2], color='k')
    plt.xlabel(col1)
    plt.ylabel(col2)
    
    plt.colorbar(contour_set_full)
    
    plt.clabel(contour_set_bands, fontsize='large', fmt='%+.0f')
    #plt.show()
    plt.savefig(col1+"To"+col2+".pdf")
