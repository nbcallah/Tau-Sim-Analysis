#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt

jtonev = 6.2415091e27
mass_n = 1.674927471e-27
grav = 9.80665
minz = -1.464413669130002
    
def plotPos(data, holdT):
    tStartDip = 350 + holdT + 20
    tEndDip = 350 + holdT + 20 + 20
    tEndInt = 350 + holdT + 20 + 100
    cutNorm = (data['energy'] > 5.875/jtonev) & (data['time'] < data['deathTime']) & (data['time'] >= tStartDip) & (data['time'] <= tEndDip) & (data['energy'] > 35/jtonev)
    weightsNorm = ((data[cutNorm]['energy']*jtonev)**1.216666667)*(np.cos(data[cutNorm]['theta'])**0.25)
    zetas = (1.0 - 0.5/(1 + np.exp(-1000*data[cutNorm]['x']))
             - np.sqrt(
                 data[cutNorm]['x']**2
                 +
                 (
                     np.sqrt(data[cutNorm]['y']**2 + (data[cutNorm]['z']-data[cutNorm]['zoff'])**2)
                     - 0.5 -0.5/(1 + np.exp(-1000*data[cutNorm]['x'])))
                 **2)
            )
    print(np.mean(zetas))
    plt.clf()
    plt.hist(zetas, bins=100, weights=weightsNorm)
    plt.show()
    
def plotSpect(data, holdT):
    tStartDip = holdT + 20
    tEndDip = holdT + 20 + 20
    tEndInt = holdT + 20 + 120
    cutNorm = (data['time'] < data['deathTime']) & (data['time'] >= tStartDip) & (data['time'] <= tEndInt)
    weightsNorm = ((data[cutNorm]['energy']*jtonev)**1.216666667)*(np.cos(data[cutNorm]['theta'])**0.25)
    plt.clf()
    plt.hist(data[cutNorm]['energy']*jtonev, bins=100, weights=weightsNorm)
    plt.show()
    plt.clf()
    plt.hist(data[cutNorm]['theta'], bins=100, weights=weightsNorm)
    plt.show()
    plt.clf()
    plt.hist(data['energy']*jtonev, bins=100)
    plt.show()
    plt.clf()
    plt.hist(data['theta'], bins=100)
    plt.show()

def plotEDiff(data, holdT):
    tStartDip = holdT + 20
    tEndDip = holdT + 20 + 20
    tEndInt = holdT + 20 + 120
    cutNorm = (data['time'] < data['deathTime']) & (data['time'] >= tStartDip) & (data['time'] <= tEndInt)
    plt.clf()
    plt.hist((data[cutNorm]['energy'] - data[cutNorm]['eStart'])*jtonev, bins=100)
    plt.yscale('log')
    plt.show()

def countUCN(data, holdT):
    tStart = holdT + 20
    tEnd = holdT + 20 + 120
    cut = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    return np.sum(cut)

def avgT(data, holdT):
    tStart = holdT + 20
    tEnd = holdT + 20 + 120
    cut = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    return np.mean(data[cut]['time'])

def plotDip2(data, holdT):
    tStart = holdT + 20
    tEndDip2 = holdT + 20 + 20
    tEnd = holdT + 20 + 120
    cut = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    cutDip2 = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEndDip2)
    plt.clf()
    plt.hist(data[cut]['time'], bins=np.arange(tStart, tEnd+1, 1))
    plt.hist(data[cutDip2]['time'], bins=np.arange(tStart, tEnd+1, 1), color='red')
    plt.show()

def ratDip2(data, holdT):
    tStart = holdT + 20
    tEndDip2 = holdT + 20 + 20
    tEnd = holdT + 20 + 120
    cut = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    cutDip2 = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEndDip2)
    sumDip1 = np.sum(cutDip2)
    norm = np.sum(cut)
    return (sumDip1/norm, sumDip1/norm*np.sqrt(1/sumDip1 + 1/norm))

def sumCleanChk(data, holdT):
    tStart = holdT
    tEnd = holdT + 20
    cut = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    return np.sum(cut)

def plotClean(data, holdT):
    cut = (data['zoff'] == -2) & (data['time'] < 0)
    plt.clf()
#    plt.hist(data[cut]['time'] + data[cut]['settlingTime'], bins=np.arange(0,50,1))
#    plt.hist(data[cut]['time'] + 50, bins=np.arange(0,100,1))
#    plt.yscale('log')
#    plt.hist(data[cut]['energy']/(mass_n*grav)+(minz+1.5), bins=100)
#    plt.hist2d(data[cut]['energy']/(mass_n*grav)+(minz+1.5), bins=100)
#    plt.show()

def countHeated(data, holdT):
    eClean = mass_n*grav*(0.38-(minz+1.5))
    cut = (data['zoff'] > 0) & (data['eStart'] < eClean) & (data['energy'] > eClean)
    return(np.sum(cut))

def countHeatedLost(data, holdT):
    eClean = mass_n*grav*(0.38-(minz+1.5))
    cut = (data['zoff'] < -2.5) & (data['eStart'] < eClean)
    return(np.sum(cut))

def countHeatedLostCleaner(data, holdT):
    eClean = mass_n*grav*(0.38-(minz+1.5))
    cut = (data['zoff'] == -3) & (data['eStart'] < eClean) & (data['energy'] > eClean)
    return(np.sum(cut))

def countHeatedLostEscape(data, holdT):
    eClean = mass_n*grav*(0.38-(minz+1.5))
    cut = (data['zoff'] == -4) & (data['eStart'] < eClean) & (data['energy'] > eClean)
    return(np.sum(cut))

def countHeatedLostHigh(data, holdT):
    eClean = mass_n*grav*(0.38+0.05-(minz+1.5))
    cut = (data['zoff'] == -3)
    plt.clf()
#    plt.hist(data[cut]['energy']*jtonev, bins=100)
#    plt.hist(data[cut]['z'], bins=100)
#    plt.hist(data[cut]['time'], bins=100)
#    plt.hist2d(data[cut]['x'], data[cut]['y'], bins=50)
    plt.hist(data[cut]['energy']/(mass_n*grav)+(minz+1.5), bins=100)
    plt.show()
    return(np.sum(cut))

def countUncleanedLost(data, holdT):
    eCleanHigh = mass_n*grav*(0.38+0.05-(minz+1.5))
    cut = (data['zoff'] == -3) & (data['eStart'] > eCleanHigh)
    return(np.sum(cut))

def countUncleaned(data, holdT):
    eClean = mass_n*grav*(0.38-(minz+1.5))
    cleanCut = (data['zoff'] == -2) & (data['eStart'] > eClean)
    bandCut = (data['eStart'] > eClean)
    return(np.sum(bandCut) - np.sum(cleanCut))

def countUncleanedHigh(data, holdT):
    eClean = mass_n*grav*(0.38+0.05-(minz+1.5))
    cleanCut = (data['zoff'] == -2) & (data['eStart'] > eClean)
    bandCut = (data['eStart'] > eClean)
    return(np.sum(bandCut) - np.sum(cleanCut))

np.random.seed(2303616184)

dt = np.dtype([('rLenFront', np.uint32, (1)),
               ('energy', np.float64, (1)),
               ('theta', np.float64, (1)),
               ('time', np.float64, (1)),
               ('settlingTime', np.float64, (1)),
               ('eperp', np.float64, (1)),
               ('x', np.float64, (1)),
               ('y', np.float64, (1)),
               ('z', np.float64, (1)),
               ('zoff', np.float64, (1)),
               ('nhit', np.int32, (1)),
               ('nhitBot', np.int32, (1)),
               ('nhitTop', np.int32, (1)),
               ('eStart', np.float64, (1)),
               ('deathTime', np.float64, (1)),
               ('rLenBack', np.uint32, (1))])

shortData = np.fromfile("../datafiles/fulltau/fullTau_Short.bin", dtype=dt)
assert(np.all(shortData['rLenFront']==shortData[0]['rLenFront']))
longData = np.fromfile("../datafiles/fulltau/fullTau_Long.bin", dtype=dt)
assert(np.all(longData['rLenFront']==longData[0]['rLenFront']))

#shortDataNoNan = shortData[~np.isnan(shortData['energy'])]
#longDataNoNan = longData[~np.isnan(longData['energy'])]

#plotEDiff(shortData, 20)
#plotEDiff(longData, 1400)

#plotClean(shortData, 20)
#plotClean(longData, 1400)

print("uncleaned   uncleanedHigh   heated   uncleanedLost   heatedLost  heatedLostCleaner  heatedLostEscape  heatedLostHigh")
print(countUncleaned(shortData, 20), countUncleanedHigh(shortData, 20), countHeated(shortData, 20), countUncleanedLost(shortData, 20), countHeatedLost(shortData, 20), countHeatedLostCleaner(shortData, 20), countHeatedLostEscape(shortData, 20), countHeatedLostHigh(shortData, 20))
print(countUncleaned(longData, 1400), countUncleanedHigh(longData, 1400), countHeated(longData, 1400), countUncleanedLost(longData, 1400), countHeatedLost(longData, 1400), countHeatedLostCleaner(longData, 1400), countHeatedLostEscape(longData, 1400), countHeatedLostHigh(longData, 1400))

numShort = countUCN(shortData, 20)
timeShort = avgT(shortData, 20)

numLong = countUCN(longData, 1400)
timeLong = avgT(longData, 1400)

#plotDip2(shortData, 20)
#plotDip2(longData, 1400)

print(ratDip2(shortData, 20), ratDip2(longData, 1400))

print(sumCleanChk(shortData, 20), sumCleanChk(longData, 1400))

lt = (timeLong-timeShort)/(np.log(numShort/numLong))
lterr = np.sqrt(1/numShort + 1/numLong)*(timeLong-timeShort)/(np.log(numShort/numLong))**2
print(timeShort, timeLong, numShort, numLong, (timeLong-timeShort), lt, lterr)