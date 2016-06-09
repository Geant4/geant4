# $Id: extractShelve.py,v 1.3 2003/07/04 13:26:31 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-09-00 $
# -------------------------------------------------------------------
#
import string
import math
import os
import shelve
import string
import detector



# some functions for calculating the FOM

def getR2(measure):
    "Error squared of a measure."
    v =  measure.variance
    m2 = math.pow(measure.mean, 2)
    n = measure.entries

    r2 = 1e+10
    if n>0 and m2 > 0:
        r2 = v/m2/n
    
    return r2

def getFOMbyName(shelveName, tallyName, bin):
    "FOM by shelve name, tally name and bin of energy region"
    she = shelve.open (shelveName, "r")
    measure = she[tallyName].measures[bin]
    R2= getR2(measure)
    T = she["runTime"]
    return 1./(T*R2)

def getFOM(aShelve, tallyName, bin):
    "FOM given a shelve, tally name and energy bin."
    measure = aShelve[tallyName].measures[bin]
    R2= getR2(measure)
    T = aShelve["runTime"]
    return 1./(T*R2)










# functions to retrive information from the shelve

def infoShelve(file):
    "Print information stored in a shelve."
    she = shelve.open(file, "r")
    print file
    for k in she.keys():
        value = she[k]
        if value.__class__.__module__ == "__builtin__" :
            print "   ", k, ":", value
            
def lsShelve(path):
    "List shelve information  found in and under the given path. "
    if os.path.isfile(path):
        if string.find(path,".shelve") > -1:
            infoShelve(path)
    else:
        if os.path.isdir(path):
            files = os.listdir(path)
            for file in files:
                f=path + "/" + file
                lsShelve(f)
                


def getFluxInRegion(she, dist="00", detType="ring", region=1):
    "Get the flux in a detector."
    tally = she["detector_" + dist + "Tally"]
    m = tally.measures[region]
    rawFlux = m.sum
    mean = m.mean
    var = m.getVariance()
    n = m.entries
    errRel = -1.
    if n>0 and var >0 and mean > 0:
        errRel = math.sqrt(var/n)/mean
    nPeakNeutrons = she["generatorTally"].measures[1].sum
    scale = detector.detScale(nPeakNeutrons, she["energy"], 0)
    return rawFlux * scale / detector.DetectorVolume[detType][dist], errRel



def getFluxes(she, detType = "ring"):
    "Get all fluxes and FOM in a shelve."
    fluxName = "flux_" + she["energy"]+"_"+\
               she["shieldWidth"]
    keys = she.keys()
    fluxes = {}
    for k in keys:
        if string.find(k,"detector_") > -1:
            sp = string.split(k,"detector_")
            sp = string.split(sp[1],"Tally")
            dist = sp[0]
            fluxNameD = fluxName + "_" + dist
            for i in range(1,3):
                fom = getFOM(she, k, i)
                f = fluxNameD + "_" + "%(i)d" % vars()
                flux, errRel = getFluxInRegion(she, dist, detType, i)
                errRel*=100
                sFlux = '%1.2E' % (flux)
                sErrRel = '%1.0f' % (errRel)
                if errRel < 10.0:
                    sErrRel = '%1.1f' % (errRel)
                print f, sFlux, ", errRel:", sErrRel, "   FOM:", fom
#                fluxes[f] =  sFlux
    return fluxes
