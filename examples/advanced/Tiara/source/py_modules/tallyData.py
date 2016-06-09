# $Id: tallyData.py,v 1.2 2003/06/16 17:06:44 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-08-00 $
# -------------------------------------------------------------------
#

class MeasureData:
    def __init__(self, entries, mean, sum, sumSquared, variance):
        self.entries = entries
        self.mean = mean
        self.sum = sum
        self.sumSquared = sumSquared
        self.variance = variance

    def __iadd__(self,m):
        self.entries += m.entries
        self.sum += m.sum
        self.sumSquared += m.sumSquared
        self.mean = self.getMean()
        self.variance = self.getVariance()
        return self
    
    def __add__(self,b):
        c = MeasureData(self.entries,
                        self.mean,
                        self.sum,
                        self.sumSquared,
                        self.variance)
        c+=b
        return c

    def getMean(self):
        return 1.0 * self.sum / self.entries

    def getVariance(self):
        n = 0
        f = 0
        if self.entries > 1:
            mean = self.getMean()
            n = 1.0 * self.entries/(self.entries -1)
            f = 1.0 * self.sumSquared/self.entries - mean*mean
        return n * f
        
    

class TallyData:
    def __init__(self, binEdges, measures):
        self.binEdges = binEdges
        self.measures = measures

    def addMeasures(self, tally):
        if tally.binEdges != self.binEdges:
            print "TallyData.addMeasures: tally.binEdges != self.binEdges"
        else:
            for i in range(len(self.measures)):
                self.measures[i] += tally.measures[i]
                
        



def createTallyDat(tally):
    measures = []
    for i in range(tally.size()):
        m = tally.measure(i)
        measures.append(MeasureData(m.GetEntries(),
                                    m.GetMean(),
                                    m.GetSum(),
                                    m.GetSumSquared(),
                                    m.GetVariance()))
    t = TallyData(tally.binEdges(), measures)
    return t
