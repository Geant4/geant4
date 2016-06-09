# $Id: myUtils.py,v 1.4 2003/06/20 12:41:06 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-06-00 $
# -------------------------------------------------------------------
#
import CLHEP
import G4Kernel

import string
import os
import time
import shelve

G4analysisUse = os.environ.has_key("G4ANALYSIS_USE")
if G4analysisUse:
    import myLiz


def createParallelSampler(impGeo, impScorer):
    parallelSampler = G4Kernel. \
                      G4ParallelGeometrySampler(\
            impGeo.getWorldVolume(), "neutron")
    parallelSampler.PrepareScoring(impScorer)
    if impGeo.base > 1:
        istore = impGeo.getImportanceStore()
        parallelSampler.PrepareImportanceSampling(istore,None)

    parallelSampler.Configure()
    return parallelSampler



def printImpTable(impScorer, istore = ""):
    iScMap = impScorer.GetMapGeometryCellCellScorer()
    table = ""
    if istore:
        table = G4Kernel.G4ScoreTable(istore)
    else:
        table = G4Kernel.G4ScoreTable()
    table.Print(iScMap)



def saveResults(tApp, path, shelveName, impScorer, impGeo):
    myShelve = shelve.open(path + "/" + shelveName)
    printImpTable(impScorer, impGeo.getImportanceStore())
    table = G4Kernel.G4ScoreTable()
    table.Print(tApp.cellScorerStore.GetMapGeometryCellCellScorer())
    if G4analysisUse:
        storePath = path + "/" + myShelve["xmlStoreName"]
        persistantStore = myLiz.tf.create(storePath,"xml",0,1)
        coppyTrees(tApp.tree, persistantStore)
        persistantStore.commit()
        persistantStore.close()
        print "wrote to store: ", storePath

    tApp.fillShelve(myShelve)
    myShelve.close()
    print "wrote shelve: ", shelveName


def coppyTrees(srcTree, cpTree):
    mntPoint = "cpTreeMountPoint"
    srcTree.mkdir("/" + mntPoint)
    srcTree.mount("/" + mntPoint, cpTree, "/")
    objectNames = srcTree.listObjectNames()
    objectTypes = srcTree.listObjectTypes()
    for i in range(len(objectNames)):
        name = objectNames[i]
        type = objectTypes[i]
        if string.find(name, mntPoint) < 0:
            srcTree.cp (name,"/" + mntPoint +"/")

min = 6000
hour = 60 * min
day = 24 * hour


def getTotalTime(totalTime):
    days = totalTime / day
    dayRest = totalTime - days*day
    hours = dayRest / hour
    minrest = dayRest - hours*hour
    minutes = minrest / min

    tString = ""
    if days > 0:
        tString += "%(days)dd" % vars()
    if hours > 0:
        tString += "%(hours)dh" % vars()
    if minutes > 0:
        tString += "%(minutes)dm" % vars()

    return tString


def getStoreName():
    
    Y, M, D, h, m, s, wd, jd, ds  = time.localtime()
    hostName = os.environ["HOST"]

    storeName = "tiara-" + \
                "%(Y)d" % vars() + "_" + \
                "%(M)d" % vars() + "_" + \
                "%(D)d" % vars() + "_" + \
                "%(h)d" % vars() + "_" + \
                "%(m)d" % vars() + "_" + \
                "%(s)d" % vars() + "_" + \
                hostName 
    
    return  storeName


def getConfigurationInfo(impGeo, experiment, physicsList, totalTime,
                         comment):
    configInfo = {}


    configInfo["energy"] = experiment.energy

    configInfo["shieldMaterial"] = experiment.shieldMaterial


    width = experiment.shieldWidth / CLHEP.cm
    s_width = "%(width)d" % vars()
    configInfo["shieldWidth"] = s_width 
    

    base = impGeo.base
    simp = "no"
    if base > 1:
        simp = "ImpBase_%(base)f" % vars()


    configInfo["biasing"] = simp


    configInfo["physListName"] = physicsList.getName()


    particles = ""

    configInfo["minEnergyCut"] = experiment.particleCut
        

    stime = "no"
    if totalTime > 0:
        stime = getTotalTime(totalTime)


    configInfo["timeLimit"] = stime


    configInfo["impGeoName"] = impGeo.nameExt

 
    configInfo["hostName"] = os.environ["HOST"]



    configInfo["comment"] = comment


    return  configInfo




def setConfigInfo(tree, confInfo, af):
    hf = af.createHistogramFactory(tree)
    h = hf.createHistogram1D("configInfo","configInfo",1,0,1)
    anno = h.annotation()
    for k in confInfo:
        anno.addItem(k, confInfo[k])

    
def getConfigInfoFromTree(tree):
    configInfo = {}
    h = tree.findH1D("configInfo")
    anno = h.annotation ()
    for i in range(8, anno.size()):
        configInfo[anno.key(i)] = anno.value(i)

    return configInfo

if G4analysisUse:
    def addToXML(mergedXMLStore, xmlStore):
        objNames = xmlStore.listObjectNames()
        objTypes = xmlStore.listObjectTypes()
        hf = myLiz.af.createHistogramFactory(mergedXMLStore)
        for i in range(len(objNames)):
            oName = objNames[i]
            oType = objTypes[i]
            if oType == "IHistogram1D":
                h1 = mergedXMLStore.findH1D(oName)
                h2 = xmlStore.findH1D(oName)
                hf.add(oName,h1,h2)
        
if G4analysisUse:
    def comparableShelves(she, mergedXMLStore):
        comp = 1
        if mergedShelve["energy"] != she["energy"]:
            comp = 0
        if mergedShelve["shieldWidth"] != she["shieldWidth"]:
            comp = 0
    

def setUpMergedShelve(she, mergedShelve, mergedXMLname, shelveNameToBeAdded):
    mergedShelve["xmlStoreName"] = mergedXMLname
    mergedShelve["energy"] = she["energy"]
    mergedShelve["shieldWidth"] = she["shieldWidth"]
    mergedShelve["mergedFiles"] = shelveNameToBeAdded
    mergedShelve["runTime"] = she["runTime"]
    mergedShelve["generatorTally"] = she["generatorTally"]
    mergedShelve["sourceDetectorTally"] = she["sourceDetectorTally"]
    mergedShelve["detector_00Tally"] = she["detector_00Tally"]
    mergedShelve["detector_20Tally"] = she["detector_20Tally"]
    mergedShelve["detector_40Tally"] = she["detector_40Tally"]
    

def addToShelve(mergedShelve, she, shelveNameToBeAdded):
    mergedShelve["mergedFiles"] += shelveNameToBeAdded + " "
    mergedShelve["runTime"] += she["runTime"]
    mergedShelve["generatorTally"].addMeasures(she["generatorTally"])
    mergedShelve["sourceDetectorTally"].addMeasures(she["sourceDetectorTally"])
    mergedShelve["detector_00Tally"].addMeasures(she["detector_00Tally"])
    mergedShelve["detector_20Tally"].addMeasures(she["detector_20Tally"])
    mergedShelve["detector_40Tally"].addMeasures(she["detector_40Tally"])

if G4analysisUse:
    def mergeData(mergedName, shelveList):
        mergedShelveName = mergedName + ".shelve"
        print "createing merged shelve: ", mergedShelveName
        mergedShelve = shelve.open(mergedShelveName)
        mergedXMLname = mergedName + ".xml"
        print "createing merged XML store: ", mergedXMLname
        mergedXMLStore = myLiz.tf.create(mergedXMLname, "xml", 0, 1)

    
        for i in range(len(shelveList)):
            shelveName = shelveList[i]
            print "adding data from: ", shelveName
            she = shelve.open(shelveName,"r")
            xmlStoreName = she["xmlStoreName"]
            print "opening: ", xmlStoreName
            xmlStore = myLiz.tf.create(xmlStoreName,"xml",1,0)
            if i == 0:
                if mergedXMLStore.listObjectNames() == ():
                    setUpMergedShelve(she, mergedShelve, mergedXMLname, shelveName)
                    coppyTrees(xmlStore, mergedXMLStore)
                else:
                    if comparableShelves(she, mergedXMLStore):
                        addToShelve(mergedShelve, she, shelveName)
                        addToXML(mergedXMLStore, xmlStore)
            else:
                addToShelve(mergedShelve, she, shelveName)
                addToXML(mergedXMLStore, xmlStore)

        mergedShelve.close()
        mergedXMLStore.commit()
        mergedXMLStore.close()


def rmPath(she):
    xmlFile = she["xmlStoreName"]
    n = string.rfind(xmlFile,"/")
    if n > -1:
        xmlFile = xmlFile[n+1:]
    she["xmlStoreName"] = xmlFile

    

