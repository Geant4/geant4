# $Id: runSequence.py,v 1.4 2004/06/09 15:04:36 daquinog Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-07-00-cand-01 $
# -------------------------------------------------------------------
#
import Tiara
import myUtils
import shelve
import os
import string

class RunConfig(object):
    def __init__(self):
        self.basePath = ""
        self.tApp = ""
        self.tiaraSpecs = ""
        self.impGeo = ""
        self.impScorer = ""
        self.totalTime = ""
        self.comment = ""
    def getConfInfo(self):
        return myUtils.\
               getConfigurationInfo(self.impGeo,
                                    self.tiaraSpecs.experiment,
                                    self.tApp.physicsList,
                                    self.totalTime,
                                    self.comment)



      
class RunSequence(object):
    # methods to be used public
    def  __init__(self, runConfig, usePI = True):
        self.rc = runConfig
        self.usePI = usePI
        self.runNum = -1
        self.confInfo = self.rc.getConfInfo()
        self.storeName = myUtils.getStoreName()
        self.xmlStore = ""
        self.pathXMLName = ""
        self.shelveName = ""
        self.pathShelveName = ""
        self.randomNumberFileName = ""
        self.path = self.mkPath()

    def runNevents(self, events):
        self.runNum += 1
        self.mkNames()
        self.mkShelve()
        self.report()
        self.rc.tApp.tiaraSim.BeamOn (events)
        Tiara.saveRandomStatus(self.randomNumberFileName)
        myUtils.saveResults(self.rc.tApp,
                            self.path,
                            self.shelveName,
                            self.rc.impScorer,
                            self.rc.impGeo)

    def runLoop(self):
        while ( not \
                (self.rc.tApp.eventAction.\
                 GetTotalProcessedTime() > self.rc.totalTime)):
            self.runNevents(10000000) # dummy num. events







    # methods used privately

    def report(self):
        print "\n\nRunSequence.report:"
        print self.confInfo
        if self.usePI:
            print "the xml store will be named: "
            print " ", self.pathXMLName
        print "the shelve name: "
        print " ", self.pathShelveName
        print "\n\n"

        

    def mkPath(self):
        if not os.path.exists(self.rc.basePath):
            os.mkdir(self.rc.basePath)
        path = self.rc.basePath + "/" + self.storeName
        if not os.path.exists(path):
            os.mkdir(path)
        return path

    def mkNames(self):
        rn = self.runNum
        rns = "%(rn)d" % vars()
        rns = string.rjust(rns, 5)
        rns = string.replace(rns,' ','0')
        rId = "_run" + rns 
        if self.usePI:
            self.xmlStore = self.storeName + rId + ".xml"
            self.pathXMLName = self.path + "/" + self.xmlStore

        self.shelveName = self.storeName + rId + ".shelve"
        self.pathShelveName = self.path + "/" + self.shelveName
        self.randomNumberFileName = self.path + "/randomNumberFile" + rId
       

    def mkShelve(self):
        myShelve = shelve.open(self.pathShelveName)
        if self.usePI:
            myShelve["xmlStoreName"] = self.xmlStore
        for info in self.confInfo:
            myShelve[info] = self.confInfo[info]
        myShelve.close()

            
# end of RunSeuence        

