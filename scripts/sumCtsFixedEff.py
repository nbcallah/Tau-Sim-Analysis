#!/usr/local/bin/python3

import sys
import numpy as np

jtonev = 6.2415091e27
np.random.seed(2303616184)

if(len(sys.argv) != 6):
	sys.exit("Error! Usage: ./plotArrivalTime fname dipStartT dipStopT integralStartT integralStopT")

fname = sys.argv[1]
dipStartT = float(sys.argv[2])
dipStopT = float(sys.argv[3])
integralStartT = float(sys.argv[4])
integralStopT = float(sys.argv[5])

dt = np.dtype([('rLenFront', np.uint32, (1)),
               ('energy', np.float64, (1)),
               ('time', np.float64, (1)),
               ('eperp', np.float64, (1)),
               ('x', np.float64, (1)),
               ('y', np.float64, (1)),
               ('z', np.float64, (1)),
               ('zoff', np.float64, (1)),
               ('nhit', np.int32, (1)),
               ('nhitBot', np.int32, (1)),
               ('nhitTop', np.int32, (1)),
               ('rLenBack', np.uint32, (1))])

data = np.fromfile(fname, dtype=dt)

integral = 0
norm = 0

for row in data:
    if row['time'] >= dipStartT and row['time'] < dipStopT:
        norm += 1
    if row['time'] >= integralStartT and row['time'] < integralStopT:
        integral += 1

print(integral/norm, np.sqrt(integral)/norm)
