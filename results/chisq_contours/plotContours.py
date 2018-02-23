#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
import argparse
import sys

parser = argparse.ArgumentParser(description='Plot Chi^2 Contours')
parser.add_argument('filename', metavar='File', type=str, help='File containing chi^2 data')
parser.add_argument('cols', metavar='cols', type=str, nargs='+', help = 'Name of Columns')
#parser.add_argument('col2', metavar='col2', type=str, help='Name of Column 2')

args = parser.parse_args()

dtList = [(c, np.float64) for c in args.cols]
dtList.append(('Chi^2', np.float64))
dt = np.dtype(dtList)

dt = np.dtype([(args.col1, np.float64),
               (args.col2, np.float64),
               ('Chi^2', np.float64)])

data = np.loadtxt(args.filename, dt)

minChisqBasin = data[data['Chi^2'] < 20]

minCol1 = np.amin(minChisqBasin[args.col1])
maxCol1 = np.amax(minChisqBasin[args.col1])
minCol2 = np.amin(minChisqBasin[args.col2])
maxCol2 = np.amax(minChisqBasin[args.col2])

minChisq = np.amin(minChisqBasin['Chi^2'])
maxChisq = np.amax(minChisqBasin['Chi^2'])

#interpRegion = data[data[args.col1] > minCol1 and data[args.col1] < maxCol1 and data[args.col2] > minCol2 and data[args.col2] < maxCol2]

xs = np.linspace(minCol1, maxCol1, 200)
ys = np.linspace(minCol2, maxCol2, 200)
xGrid, yGrid = np.meshgrid(xs, ys)
zsInterp = interp.griddata((minChisqBasin[args.col1], minChisqBasin[args.col2]), minChisqBasin['Chi^2'], (xGrid, yGrid), method='linear')

#print()

plt.contour(xGrid, yGrid, zsInterp, [minChisq, minChisq+1, minChisq+2], colors='k')
#plt.contour(xGrid, yGrid, zsInterp, 15)
plt.contourf(xGrid, yGrid, zsInterp, 15, norm=plt.Normalize(vmax=maxChisq, vmin=minChisq))
#plt.scatter(minChisqBasin[args.col1], minChisqBasin[args.col2], color='k')
plt.show()