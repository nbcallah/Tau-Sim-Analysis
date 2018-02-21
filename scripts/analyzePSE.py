#!/usr/local/bin/python3

import sys
import numpy as np

def sumCts(data, dipStartT, dipStopT, integralStartT, integralStopT):
    integral = 0
    norm = 0
    for row in data:
        if row['time'] >= dipStartT and row['time'] < dipStopT:
            norm += 1
        if row['time'] >= integralStartT and row['time'] < integralStopT:
            integral += 1
    return (integral/norm, np.sqrt(integral)/norm)

jtonev = 6.2415091e27
np.random.seed(2303616184)

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

holdTs = [20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
fNames = {}
for t in holdTs:
    fNames[t] = "../datafiles/pse/pse_512k_"+str(t)+".bin"

data = np.fromfile(fNames[20], dtype=dt)
(shortCts, shortErr) = sumCts(data, 200+20, 200+20+120, 200+20+20, 200+20+40)

for t in holdTs:
    data = np.fromfile(fNames[t], dtype=dt)
    (cts, err) = sumCts(data, 200+t, 200+t+120, 200+t+20, 200+t+40)
#    print(t, cts, 0.0, err)
    print(t, cts/shortCts, 0.0, err/shortCts)
