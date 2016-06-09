# $Id: dataAcess.py,v 1.2 2003/06/16 17:06:44 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-05-02 $
# -------------------------------------------------------------------
#
import os
import shelve
import myLiz
import dpsManip
import detector



class ExperimentalData(object):
    def __init__(self):
        tiara_dir = os.environ["TIARA_BASE"]
        if not tiara_dir:
            print "dataAcess.ExperimentalData: TIARA_BASE not defined run tiara...sh first"
        dataFile = tiara_dir + "/data/expDataConverted/TiaraData2.xml"
        print "display.ExperimentalData.dataFile:" ,dataFile
        self.dataTree = myLiz.tf.create (dataFile,"xml",1,0)

    def getDataDPS(self, she, detector):
        energy  = she["energy"]
        shieldWidth = she["shieldWidth"]
        dataName = "Tiara-" + energy + "c" + shieldWidth + \
                   "-" + detector
        if shieldWidth == "25" or shieldWidth == "50":
            dataName += "a"
        dataName += ".pnt"
        print dataName
        pData = self.dataTree.findDataPointSet(dataName)
        dpsManip.setDPSErrorsToZero(pData, 0)
        return pData

    def getScaledDataaDPS(self, she, detname, df):
        pDataO = self.getDataDPS(she, detname)
        pData = dpsManip.createScaledDPS(0,
                                         pDataO,
                                         df,
                                         "scaled_" + pDataO.title (),
                                         0.000001)
        return pData



class MC_Data(object):
    def __init__(self, she, mcTree):
        self.she = she
        self.energy = self.she["energy"]
        self.shield = self.she["shieldWidth"]
        self.mcTree = mcTree
        self.coli = 0
        if self.shield == "25"  or  \
           self.shield == "50":
            self.coli = 1
        self.baseName = self.energy + "c" + self.shield +\
                        "_detector_"
        
     
    def getGeneratedHisto(self):
        name = "source_detector"
        hGen = self.mcTree.findH1D(name)
        return hGen
    
    def getMcPlot(self, detector, histo):
        mcName = "detector_" + detector + histo
        print mcName
        hMc = self.mcTree.findH1D(mcName)
        return hMc

    def getScale(self, atColiExit):
        coli = 0
        if atColiExit==1:
            coli = 0
        else:
            coli = self.coli

        nPeakNeutrons = self.she["generatorTally"].measures[1].sum
        scale = detector.detScale(nPeakNeutrons, self.energy, coli)
        print "CompPlot.getScale: scaling with:", scale
        return scale


    def getScaledMcDPS(self, atColiExit, det, df, histo = ""):
        dname = self.baseName + det.name
        if histo:
            dname += "_" + histo
    
        scale = self.getScale (atColiExit)
        hMc = self.getMcPlot (det.name, histo )
        binEdges = dpsManip.getBinEdges(hMc)
        dScaled = dpsManip.getScaledDPS(hMc, scale/det.volume,
                                        df, dname + "scaled")
        dLethScaled = dpsManip.dLogWeightDPS (dScaled,
                                              df,
                                              dname + "df_dlgE",
                                              binEdges)
        return dLethScaled


    def getScaledGeneratedDPS(self, atColiExit, det, df, histo = ""):
        dname = self.baseName + det.name
        if histo:
            dname += "_" + histo
        sourceDetectorVolume = 9.33
        scale = self.getScale (atColiExit)
        hGen = self.getGeneratedHisto()
        binEdges = dpsManip.getBinEdges(hGen)
        dGenScaled = dpsManip.getScaledDPS(hGen,
                                           scale/sourceDetectorVolume,
                                           df,
                                           dname + "GenScaled")
        dGenLethScaled = dpsManip.dLogWeightDPS(dGenScaled, df,
                                                dname + "gen, df_dlgE",
                                                binEdges)

        return dGenLethScaled



class ExpMcPlot(object):
    """Prepare and hold source information for a plot.

    Hold experimental and Monte Carlo data to plot in one diagram. 
    """
    def __init__(self, shelveName, dist, detType="ring", histo = ""):
        path = os.path.dirname(shelveName)
        shelveFile = os.path.basename(shelveName)
        self.tt = myLiz.tf.create ()
        self.df = myLiz.af.createDataPointSetFactory (self.tt)

        self.she = shelve.open(shelveName,"r")
        xmlFile = self.she["xmlStoreName"]
        if path:
            xmlFile = path + "/" + xmlFile
        self.mcTree = myLiz.tf.create (xmlFile, "xml", 1, 0)
        self.det = detector.Detector(dist,detType)
        self.expData = ExperimentalData()
        self.mcData = MC_Data(self.she, self.mcTree)

        self.pDataDPS = self.expData.getScaledDataaDPS(self.she,
                                                       self.det.name,
                                                       self.df)
        
        self.pMcDPS = self.mcData.getScaledMcDPS(atColiExit=1,
                                                 det=self.det,
                                                 df=self.df,
                                                 histo="")
        
            
        self.pGenDPS = self.mcData.getScaledGeneratedDPS(atColiExit=1,
                                                      det=self.det,
                                                      df=self.df,
                                                      histo="")

        self.regions = []


    def display(self):
        if "pl" not in dir(myLiz):
            myLiz.pf = myLiz.af.createPlotterFactory ()
            myLiz.pl = myLiz.pf.create()
        region = myLiz.pl.currentRegion()
        self.regions.append(region)
        region.plot (self.pDataDPS,"markers overlay")
        region.plot (self.pMcDPS,"markers overlay")
        myLiz.pl.refresh ()







