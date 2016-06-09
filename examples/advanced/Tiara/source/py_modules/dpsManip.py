# $Id: dpsManip.py,v 1.2 2003/06/16 17:06:44 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-07-01 $
# -------------------------------------------------------------------
#
import math

# functions to manipulate a DPS (data point set)

def getScaledDPS(hSrc, scale, df, name):
    "Create a scaled DPS from a histogramm."
    hSrc.scale(scale)
    dp = df.create(name, hSrc)
    setDPSErrorsToZero(dp, 0)
    hSrc.scale(1./scale)

    return dp


def setDPSErrorsToZero(dps, coNum):
    "Set errors of a DPS to zero."
    for i in range(dps.size()):
        p = dps.point(i)
        co = p.coordinate(coNum)
        co.setErrorMinus(coNum)
        co.setErrorPlus(coNum)

def copyDPS(dSrc, df, name):
    "Make a copy of a DPS using the data point set factory df."
    dp = df.create(name, dSrc.dimension ())
    for i in range( dSrc.size() ):
        dp.addPoint()
        pNew = dp.point(i)
        pOld = dSrc.point(i)
        for d in range(dSrc.dimension ()):
            cOld = pOld.coordinate(d)
            cNew = pNew.coordinate(d)
            cNew.setValue(cOld.value())
            em = cOld.errorMinus()
            ep = cOld.errorPlus()
            cNew.setErrorMinus(em)
            cNew.setErrorPlus(ep)
    return dp

def dLogWeightDPS(dSrc, df, name, binEdges):
    """Create a dps with the coordinate 1 scaled by 1/dlog(coordinate 0).
    """
    dp = copyDPS(dSrc, df, name)
    for i in range(len(binEdges) - 1):
        s = math.log(binEdges[i+1]) - math.log(binEdges[i])
        c = dp.point(i).coordinate (1)
        value = dSrc.point(i).coordinate(1).value()
        errorMinus = dSrc.point(i).coordinate(1).errorMinus()
        errorPlus = dSrc.point(i).coordinate(1).errorPlus()
        c.setValue(value / s)
        c.setErrorMinus(errorMinus / s)
        c.setErrorPlus(errorPlus / s)
    return dp


def createScaledDPS(coNum, pDataO, df, name, scale):
    """Create a scaled DPS from a DPS.

    The coordinate \'coNum\' of the source DPS \'pDataO\' is scaled by
    \'scale\'. The result is returned in the DPS created by  the given
    data point set factory \'df\' named \'name\'.
    """
    dp = copyDPS(pDataO, df, name)
    for i in range( dp.size() ):
        p = dp.point(i)
        co = p.coordinate(coNum)
        co.setValue(co.value()*scale)
        co.setErrorMinus(co.errorMinus()*scale)
        co.setErrorPlus(co.errorPlus()*scale)
    return dp



def getBinEdges(h):
    "Get bin edges of a histogramm."
    binEdges = []
    a = h.axis ()
    for i in range(a.bins ()):
        binEdges.append(a.binLowerEdge(i))
    binEdges.append(a.binUpperEdge(a.bins () - 1))
    return binEdges




