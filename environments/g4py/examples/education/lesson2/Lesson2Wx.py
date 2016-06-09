#!/usr/bin/python

#### 2006 Sep 26, the first draft version
## Wired not yet. g4pipe control and Popen

from Geant4 import *
import g4py.Qmaterials, g4py.NISTmaterials
import g4py.ExN03geom
import g4py.ExN03pl
import g4py.ParticleGun, g4py.MedicalBeam
import sys
from time import *
from subprocess import *
import os

# ==================================================================
# main
# ==================================================================
# ------------------------------------------------------------------
# randum number
# ------------------------------------------------------------------
rand_engine= Ranlux64Engine()
HepRandom.setTheEngine(rand_engine)
HepRandom.setTheSeed(20050830L)

# ------------------------------------------------------------------
# setup for materials
# ------------------------------------------------------------------


# NIST materials
#g4py.NISTmaterials.Construct()

# ------------------------------------------------------------------
# setup for geometry
# ------------------------------------------------------------------
# normal way for constructing user geometry

exN03geom= g4py.ExN03geom.ExN03DetectorConstruction()
gRunManager.SetUserInitialization(exN03geom)

# 2nd way, short-cut way

#g4py.ExN01geom.Construct()
#g4py.ExN03geom.Construct()

# magnetic field
#exN03geom.SetMagField(0.1 * tesla)

# ------------------------------------------------------------------
# setup for physics list
# ------------------------------------------------------------------
# normal way for constructing user physics list
exN03PL= g4py.ExN03pl.PhysicsList()
gRunManager.SetUserInitialization(exN03PL)

# 2nd way, short-cut way
#g4py.ExN01pl.Construct()


# ------------------------------------------------------------------
# setup for primary generator action
# ------------------------------------------------------------------
# normal way for constructing user physics list
#pgPGA= g4py.ParticleGun.ParticleGunAction()
#gRunManager.SetUserAction(pgPGA)
#pg= pgPGA.GetParticleGun()

# 2nd way, short-cut way
pg= g4py.ParticleGun.Construct()

# set parameters of particle gun
pg.SetParticleByName("e-")
pg.SetParticleEnergy(50.*MeV)
pg.SetParticlePosition(G4ThreeVector(-40.,0.,0.)*cm)
pg.SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.))

# medical beam
#beam= MedicalBeam.Construct()

# ------------------------------------------------------------------
# go...
# ------------------------------------------------------------------
gRunManager.Initialize()


# visualization
# OGLSX, VRML and HEPREP sceneHandlers are all created with names
gApplyUICommand("/vis/sceneHandler/create OGLSX OGLSX")
gApplyUICommand("/vis/sceneHandler/create VRML2FILE VRML")
gApplyUICommand("/vis/sceneHandler/create HepRepFile HEPREP")

#  OGLSX is the default so, viewer is created and volume is drawn
gApplyUICommand("/vis/viewer/create OGLSX oglsxviewer")
gApplyUICommand("/vis/drawVolume")
gApplyUICommand("/vis/scene/add/trajectories")

gApplyUICommand("/tracking/storeTrajectory 1")
gApplyUICommand("/vis/scene/endOfEventAction accumulate")
gApplyUICommand("/vis/scene/endOfRunAction accumulate")
gApplyUICommand("/vis/viewer/select  oglsxviewer")

# viewers VRML and Wired are tested by their envs vars
# if their envs var are set, then viewers are created and drawVolume

global heprepViewer, heprepDir, heprepName, vrmlViewer
heprepViewer = os.environ.get("G4HEPREPFILE_VIEWER")
heprepDir = os.environ.get("G4HEPREPFILE_DIR")
heprepName = os.environ.get("G4HEPREPFILE_NAME")
if heprepViewer is not None:
  gApplyUICommand("/vis/viewer/create HEPREP wired")
  gApplyUICommand("/vis/drawVolume")

# VRML viewers name is user defined
vrmlDir = os.environ.get("G4VRML_DEST_DIR")
vrmlViewer = os.environ.get("G4VRMLFILE_VIEWER")

if vrmlViewer is not None:
  gApplyUICommand("/vis/viewer/create VRML vrmlviewer")
  gApplyUICommand("/vis/drawVolume")


# test############################################################################TESSSSSSST
gApplyUICommand("/vis/viewer/select oglsxviewer")
gApplyUICommand("/vis/scene/add/trajectories")

gApplyUICommand("/tracking/storeTrajectory 1")
gApplyUICommand("/vis/scene/endOfEventAction accumulate")
gApplyUICommand("/vis/scene/endOfRunAction accumulate")
#test########################################################################################

global commandDic, commandList
commandDic = {}
commandList = []

def DumpTree(atree):

  ntree= atree.GetTreeEntry()
  ncommand= atree.GetCommandEntry()
  for i in range(1, ncommand+1):
    icommand= atree.GetCommand(i)
    command = str( icommand.GetCommandPath())
    commandList.append(command)
    nparameter= icommand.GetParameterEntries()
    pguide = ""
    for j in range(0, nparameter):
      iparam= icommand.GetParameter(j)
      pguide = pguide + "Parameter: " + str(iparam.GetParameterName())+ "   Type: " + 	        str(iparam.GetParameterType()) + '\n'
    guide = str(icommand.GetTitle()) + '\n' + pguide
    commandDic[command] = guide

  for i in range(1, ntree+1):
    itree= atree.GetTree(i)
    DumpTree(itree)


root_tree= gUImanager.GetTree()
DumpTree(root_tree)




######  wxPython GUI ##########


import wx



class ComPanel(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent, -1)

    self.g4comText = wx.TextCtrl(parent, -1, "enter Command", pos=(10,10), size=(300, 30))    
    self.g4comExec = wx.Button(parent, -1, "Execute", pos=(320,10), size=(60,30))
    self.g4comExec.Bind(wx.EVT_BUTTON, self.ExecuteCommand, self.g4comExec)
#    self.sizerE = wx.BoxSizer(wx.HORIZONTAL)
#    self.sizerE.Add(self.g4comText)
#    self.sizerE.Add(self.g4comExec)
    
#    self.sizerL = wx.BoxSizer(wx.HORIZONTAL)
    self.comListBox = wx.ListBox(parent, -1, pos=(10,50), size=(300,200), choices=commandList, style=wx.LB_SINGLE)
    self.comListBox.SetSelection(1)
    self.comListBox.Bind(wx.EVT_LISTBOX, self.ShowGuide, self.comListBox)
#    self.sizerL.Add(self.comListBox)

    self.guide = wx.TextCtrl(parent, -1, "guidance", pos=(320,50), size=(300, 200), style =wx.TE_MULTILINE)
    self.guide.Bind(wx.EVT_LISTBOX, self.ShowGuide, self.guide)

#    self.sizerL.Add(self.guide)
#    self.sizer = wx.BoxSizer(wx.VERTICAL)
#    self.sizer.Add(self.sizerE)
#    self.sizer.Add(self.sizerL)
#    self.SetSizer(self.sizer)

  def ShowGuide(self, event):
    self.guide.Clear() # how to cleat the whole text before showing the next
    g4com =str(self.comListBox.GetStringSelection())
    self.guide.WriteText( commandDic[g4com])
    self.g4comText.SetValue(g4com)

  def ExecuteCommand(self, event):
    gApplyUICommand(str(self.g4comText.GetValue()))

class VisPanel(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent, -1)

    self.visZoomIn = wx.Button(parent, -1, "Zoom In", pos=(10,100), size=(80,30))
    self.visZoomOut = wx.Button(parent, -1, "Zoom out", pos=(100,100), size=(80,30))
    self.visUp = wx.Button(parent, -1, "Up", pos=(190,100), size=(80,30))
    self.visDown = wx.Button(parent, -1, "Down", pos=(270,100), size=(80,30))
    self.visLeft = wx.Button(parent, -1, "Left", pos=(360,100), size=(80,30))
    self.visRight = wx.Button(parent, -1, "Left", pos=(450,100), size=(80,30))
    
    viewerList = ["OpenGL", "VRML", "Wired"]
    self.viewer = wx.RadioBox(self, -1, "Viewer", pos=(10,10),
                             size=(210,60), choices=viewerList, majorDimension=1, style=wx.RA_SPECIFY_ROWS)
    self.viewer.Bind(wx.EVT_RADIOBOX, self.ViewerSelected, self.viewer)
    self.viewer.SetToolTip(wx.ToolTip("Select one"))
    self.viewer.SetSelection(0)
    if vrmlViewer == None: self.viewer.EnableItem(1, False)
    if heprepViewer == None: self.viewer.EnableItem(2, False)    
#    self.sizer = wx.BoxSizer(wx.HORIZONTAL)
#    self.sizer.Add((20,30))
#    self.sizer.Add(self.visZoomIn)
#    self.sizer.Add(self.visZoomOut)
#    self.sizer.Add(self.visUp)
#    self.sizer.Add(self.visDown)
#    self.sizer.Add(self.visRight)
#    self.sizer.Add(self.visLeft)
#    self.SetSizer(self.sizer)
    
#    self.visZoomIn = wx.Button(parent, -1, "Zoom In", pos=(10,70), size=(80,30))
#    self.visZoomOut = wx.Button(parent, -1, "Zoom out", pos=(100,70), size=(80,30))
#    self.visUp = wx.Button(parent, -1, "Up", pos=(10,50), size=(190,70))
#    self.visDown = wx.Button(parent, -1, "Down", pos=(100,50), size=(270,70))
#    self.visLeft = wx.Button(parent, -1, "Left", pos=(10,90), size=(360,70))
#    self.visRight = wx.Button(parent, -1, "Left", pos=(100,90), size=(450,70))

                        

    self.visZoomIn.Bind(wx.EVT_BUTTON, self.cmdExpand, self.visZoomIn)
    self.visZoomOut.Bind(wx.EVT_BUTTON, self.cmdShrink, self.visZoomOut)
    self.visUp.Bind(wx.EVT_BUTTON, self.cmdUp, self.visUp)
    self.visDown.Bind(wx.EVT_BUTTON, self.cmdDown, self.visDown)
    self.visRight.Bind(wx.EVT_BUTTON, self.cmdRight, self.visRight)
    self.visLeft.Bind(wx.EVT_BUTTON, self.cmdLeft, self.visLeft)

  def cmdExpand(self, event):
    gApplyUICommand("/vis/viewer/zoom 1.2")
  def cmdShrink(self, event):
    gApplyUICommand("/vis/viewer/zoom 0.8")
  def cmdUp(self, event):
    gApplyUICommand("/vis/viewer/pan "   + " 0.  10. mm")
  def cmdDown(self, event):
    gApplyUICommand("/vis/viewer/pan " +  " 0. -10.  mm")
  def cmdRight(self, event):
    gApplyUICommand("/vis/viewer/pan " +  " -10. 0.  mm")
  def cmdLeft(self, event):
    gApplyUICommand("/vis/viewer/pan "   + " 10. 0. mm")

  def ViewerSelected(self, event):
    self.viewerName = event.GetString()
    
    if self.viewerName == "OpenGL":
      gApplyUICommand("/vis/viewer/select oglsxviewer")
      gApplyUICommand("/vis/scene/add/trajectories")

      gApplyUICommand("/tracking/storeTrajectory 1")
      gApplyUICommand("/vis/scene/endOfEventAction accumulate")
      gApplyUICommand("/vis/scene/endOfRunAction accumulate")

    if self.viewerName == "VRML":
      gApplyUICommand("/vis/viewer/select vrmlviewer")
      gApplyUICommand("/vis/scene/add/trajectories")

      gApplyUICommand("/tracking/storeTrajectory 1")
      gApplyUICommand("/vis/scene/endOfEventAction accumulate")
      gApplyUICommand("/vis/scene/endOfRunAction accumulate")

    if self.viewerName == "Wired":
      
      gApplyUICommand("/vis/viewer/select wired")
      gApplyUICommand("/vis/scene/add/trajectories")

      gApplyUICommand("/tracking/storeTrajectory 1")
      gApplyUICommand("/vis/scene/endOfEventAction accumulate")
      gApplyUICommand("/vis/scene/endOfRunAction accumulate")

# everytime wired is chosen, a new instance of wired is created
# to reuse single wired, g4pipe.poll() must be checked BEFORE the SECOND Popen
#      if g4pipe.poll() == None:
      g4pipe=Popen(heprepViewer+ " -file " + heprepDir+"/" +heprepName +".heprep", shell=True)



# not used
class MyText(wx.StaticText):
    def __init__(self, parent, Text):
        wx.StaticText.__init__(self, parent, -1, Text, pos=(20,20))
        self.Bind(wx.EVT_LEFT_UP, self.ChangeColor)
          
    def ChangeColor(self, event):
        TheColour = self.GetForegroundColour()
        if TheColour == (0,0,0):
             self.SetLabel("Simulation is running!")
             self.SetForegroundColour("red")
        else:
             self.SetLabel("Click me to start a run")
             self.SetForegroundColour("black")

# to be used to choose materials and particles
# myList is a list of keys() of a dictionary
# f.e., materials are Python objects with their names as their keys
class SelectOne(wx.RadioBox):
    def __init__(self, parent, myTitle, myList):
        wx.RadioBox.__init__(self, parent, -1, myTitle, wx.DefaultPosition,
                             wx.DefaultSize, myList, 5, wx.RA_SPECIFY_ROWS)
#        self.Bind(wx.EVT_RADIOBOX, self.Selected)
        self.SetToolTip(wx.ToolTip("Select one"))
        self.SetSelection(0)
#    used only to test  Bind and getValue
#    def Selected(self, event):
#        self.selected = event.GetString()


# used to %3.3f floating point number
# energy and length unit is given by unitList which must be a dictionary
#           with units of Python objects and their name as their keys

class FloatCounter(wx.Panel):
    def __init__(self, parent, myTitle, unitList):
        wx.Panel.__init__(self, parent, -1)

        self.sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer.Add((10, -1))                
        self.sizer.Add(wx.StaticText(parent, -1, myTitle, wx.DefaultPosition, (100, -1)))

        self.intPart = wx.SpinCtrl(parent, -1, "", wx.DefaultPosition, (60,-1))
        self.intPart.SetRange(0,999)
        self.intPart.SetValue(1)
        self.intPart.Bind(wx.EVT_SPINCTRL, self.SetFloat)

        self.manPart = wx.SpinCtrl(parent, -1, "", wx.DefaultPosition, (50, -1))
        self.manPart.SetRange(0,999)
        self.manPart.SetValue(0)
        self.manPart.Bind(wx.EVT_SPINCTRL, self.SetFloat)

        self.unitSel = wx.Choice(parent, -1, wx.DefaultPosition, (90, -1),  unitList)
        self.unitSel.Bind(wx.EVT_CHOICE, self.SetFloat, self.unitSel)
        self.unitSel.SetSelection(2)

        self.valAndUnit = wx.TextCtrl(parent, -1, "value unset", wx.DefaultPosition, (150, -1))

        self.sizer.Add(self.valAndUnit)
        self.sizer.Add((10, -1))        
        self.sizer.Add(wx.StaticText(parent, -1, "  ", wx.DefaultPosition, (30, -1)))        
        self.sizer.Add(self.intPart)
        self.sizer.Add(wx.StaticText(parent, -1, ".", wx.DefaultPosition, (10, -1)))
        self.sizer.Add(self.manPart)
        self.sizer.Add((5,-1))
        self.sizer.Add(self.unitSel)
        self.SetSizer(self.sizer)

    def SetFloat(self, event):
        self.theValueStr = str(self.intPart.GetValue()) + "." + str(self.manPart.GetValue()) + "  "
        self.theValue = float (self.intPart.GetValue()) + float(self.manPart.GetValue())/ 1000.
        self.theUnit =  self.unitSel.GetStringSelection()
        theText = "%.3f" % (self.theValue) + " " + self.theUnit
        self.valAndUnit.SetValue(theText)
# user may edit the Entry, so finnaly this value must be got
#        but it doesn't work for length but work for energy (gApplyUIcommnd)
        theTextWithStar = "%.3f" % (self.theValue) + " * " + self.theUnit

# special class for this example to set/unset processes
class Processes(wx.Panel):
    def __init__(self, parent, myTitle, myList):
        wx.Panel.__init__(self, parent, -1)
        self.processCheck = {}
        self.sizer = wx.FlexGridSizer(rows=5)
        self.sizer.AddGrowableRow(1)
        for item in myList:
            self.processCheck[item] = wx.CheckBox(parent, -1, item)
            self.processCheck[item].SetValue(True)
#            self.processCheck[item].Bind(wx.EVT_CHECKBOX, self.CheckedProcess)
            self.sizer.Add(self.processCheck[item],0,wx.EXPAND)

        self.SetSizer(self.sizer)
        self.SetBackgroundColour('green')
        self.processState = {}
        self.myList = myList

# test only
#    def CheckedProcess(self, event):
#        self.processName = event.GetEventObject().GetLabel()
#        self.processState = event.GetEventObject().GetValue()


# slider to set an integer value
# title is shown
class Adjuster(wx.Panel):
    def __init__(self, parent, myTitle, minVal, maxVal, initVal):
        wx.Panel.__init__(self, parent, -1)
        self.sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer.Add((10, -1))
        self.sizer.Add(wx.StaticText(parent, -1, myTitle, wx.DefaultPosition, (100, -1)))
        self.slider = wx.Slider(parent, -1, initVal, minVal, maxVal,
                           wx.DefaultPosition, (300, -1),
                           wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
        self.slider.SetPageSize(1)
        self.sizer.Add(self.slider)
        self.SetSizer(self.sizer)

# test only        
#        self.Bind(wx.EVT_SLIDER, self.Adjusted)

#    def Adjusted(self, event):
#        print self.GetValue()


########################## no use now
class Counter(wx.SpinCtrl):
    def __init__(self, parent, myTitle, minVal, maxVal, initVal):
        wx.SpinCtrl.__init__(self, parent, -1, "", wx.DefaultPosition, wx.DefaultSize, wx.TE_RIGHT)
        self.SetRange(minVal, maxVal)
        self.SetValue(initVal)
        self.Bind(wx.EVT_SPINCTRL, self.Adjusted)
    def Adjusted(self, event):
        print self.GetValue()

############################

# main class to instantiate the above classes and pack them using nested sizers
g4pipe=0

class MyApp(wx.Frame):
    def __init__(self):
         wx.Frame.__init__(self, None, -1, "Geant4Py")
     	 self.nb = wx.Notebook(self, -1, wx.DefaultPosition, wx.DefaultSize,
                             style=
                             wx.NB_TOP # | wx.NB_MULTILINE
                             #wx.NB_BOTTOM
                             #wx.NB_LEFT
                             #wx.NB_RIGHT
                             )
         self.nb.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged)
    	 self.nb.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.OnPageChanging)

         panel = wx.Panel(self.nb)
         self.nb.AddPage(panel, "ExampleN03")
	 commandPanel = wx.Panel(self.nb)
	 self.nb.AddPage(commandPanel, "Geant4 Commands")
	 comP = ComPanel(commandPanel)
	 gxsizer = wx.BoxSizer(wx.HORIZONTAL)
	 gxsizer.Add(comP)
	 commandPanel.SetSizer(gxsizer)
	 

         visualizationPanel = wx.Panel(self.nb)
         self.nb.AddPage(visualizationPanel, "Vis Commands")
         visP = VisPanel(visualizationPanel)
         vxsizer = wx.BoxSizer(wx.HORIZONTAL)
         vxsizer.Add(visP)
         visualizationPanel.SetSizer(vxsizer)
	 
#  outmost sizer in the vertical direction
         bxsizer = wx.BoxSizer(wx.VERTICAL)
#   nested sizer in the horizontal direction
         bysizer = wx.BoxSizer(wx.HORIZONTAL)


         self.runStart = wx.Button(panel, -1, "  Run Start", wx.DefaultPosition, wx.DefaultSize)
         self.Bind(wx.EVT_BUTTON, self.RunStart, self.runStart)
         bxsizer.Add(self.runStart, 0, wx.ALL)
# widgets
         
         absorberMaterialList = ['Aluminium', 'Lead']
         self.theAbsorberMaterial = SelectOne(panel, "Absorber Materials", absorberMaterialList)
         gapMaterialList = ["liquidArgon","Scintillator", "Air", "Aerogel",  "Galactic"]
         self.theGapMaterial = SelectOne(panel, "Gap Materials", gapMaterialList)
  
         particleList = ["proton", "gamma", "e-",  "e+", "mu-", "mu+"]
         self.theParticle = SelectOne(panel, "Particles", particleList)

         self.processList = ["phot", "compt", "conv", "msc", "eIoni", "eBrem", "annihil","muIoni", "muBrems", "hIoni"]
         self.theProcesses = Processes(panel, "Processes", self.processList)

         self.eventNo = Adjuster(panel, "Number of Events", 1 , 100 , 1)
	 self.layerNo = Adjuster(panel, "Number of Layers", 1, 10, 10)
	 
         self.lengthUnit = {'micrometer':micrometer, 'mm':mm, 'cm':cm, 'm':m}
         self.absorberThickSpin = FloatCounter(panel, "Absorber Thickness", self.lengthUnit.keys())
         self.gapThickSpin = FloatCounter(panel, "Gap Thickness", self.lengthUnit.keys())
         self.sizeYZSpin = FloatCounter(panel, "Section Size", self.lengthUnit.keys())
         self.cutLengthSpin = FloatCounter(panel, "Cut Length", self.lengthUnit.keys())
         self.magneticUnit = {'Tesla':tesla, 'gauss':gauss, 'kilogauss':kilogauss}
         self.magneticFieldSpin = FloatCounter(panel, "Magnetic Field", self.magneticUnit.keys())

         self.energyUnit = { 'keV':keV, 'MeV':MeV, 'GeV':GeV, 'TeV':TeV, 'PeV':PeV}
         self.energySpin = FloatCounter(panel, "incident beam energy", self.energyUnit.keys())
	 

# now sizers

         bxsizer.Add(wx.StaticLine(panel), 0, wx.EXPAND|wx.TOP|wx.BOTTOM, 5)
         bysizer.Add((10, -1))
         bysizer.Add(self.theAbsorberMaterial, 0, wx.EXPAND, 10)
         bysizer.Add((10, -1))
         bysizer.Add(self.theGapMaterial, 0, wx.EXPAND, 10)
         bysizer.Add((10, -1))
         bysizer.Add(self.theParticle, 0, wx.EXPAND, 10)
         bysizer.Add((10, -1))         
         bysizer.Add(self.theProcesses.sizer, 0, wx.EXPAND, 10)
         
         bxsizer.Add(bysizer, 0, wx.EXPAND)

         bxsizer.Add(wx.StaticLine(panel), 0, wx.EXPAND|wx.TOP|wx.BOTTOM, 5)
         bxsizer.Add(self.layerNo.sizer, 0, wx.EXPAND)
         bxsizer.Add(self.absorberThickSpin.sizer, 0, wx.EXPAND)         
         bxsizer.Add(self.gapThickSpin.sizer, 0, wx.EXPAND)         
         bxsizer.Add(self.sizeYZSpin.sizer, 0, wx.EXPAND)         
         bxsizer.Add(wx.StaticLine(panel), 0, wx.EXPAND|wx.TOP|wx.BOTTOM, 5)         
	 
	 bxsizer.Add(self.energySpin.sizer, 0, wx.EXPAND)
         bxsizer.Add(wx.StaticLine(panel), 0, wx.EXPAND|wx.TOP|wx.BOTTOM, 5)
         bxsizer.Add(self.cutLengthSpin.sizer, 0, wx.EXPAND)         
         bxsizer.Add(wx.StaticLine(panel), 0, wx.EXPAND|wx.TOP|wx.BOTTOM, 5)         
         bxsizer.Add(self.magneticFieldSpin.sizer, 0, wx.EXPAND)         

         bxsizer.Add(wx.StaticLine(panel), 0, wx.EXPAND|wx.TOP|wx.BOTTOM, 5)
         bxsizer.Add(self.eventNo.sizer, 0, wx.EXPAND)

#         self.solid = EZgeom.G4EzVolume.GetSold(water_phantom)
#         gControlExecute("oglx.mac")
         
         panel.SetSizer(bxsizer)
         bxsizer.Fit(self)
         bxsizer.SetSizeHints(self)
    def OnPageChanged(self, event):
        old = event.GetOldSelection()
        new = event.GetSelection()
        sel = self.nb.GetSelection()
        event.Skip()

    def OnPageChanging(self, event):
        old = event.GetOldSelection()
        new = event.GetSelection()
        sel = self.nb.GetSelection()
        event.Skip()

    def RunStart(self, event):
	    
	absorberTh = self.absorberThickSpin.theValue * self.lengthUnit[self.absorberThickSpin.theUnit]/2.0
        gapTh = self.gapThickSpin.theValue * self.lengthUnit[self.gapThickSpin.theUnit]/2.0
        yzSize = self.sizeYZSpin.theValue * self.lengthUnit[self.sizeYZSpin.theUnit]
        print "BUG"
        cutLen = self.cutLengthSpin.theValue * self.lengthUnit[self.cutLengthSpin.theUnit]
        magF = self.magneticFieldSpin.theValue * self.magneticUnit[self.magneticFieldSpin.theUnit]
        
	exN03geom.SetNbOfLayers(self.layerNo.slider.GetValue())
      	exN03geom.SetAbsorberMaterial(str(self.theAbsorberMaterial.GetStringSelection()))
      	exN03geom.SetAbsorberThickness(absorberTh)
      	exN03geom.SetGapMaterial(str(self.theGapMaterial.GetStringSelection()))
      	exN03geom.SetGapThickness(gapTh)
      	exN03geom.SetCalorSizeYZ(yzSize)
      	position = -self.layerNo.slider.GetValue() * ( absorberTh + gapTh )*1.2

      	exN03geom.UpdateGeometry()
      	exN03PL.SetDefaultCutValue(cutLen)
      	exN03PL.SetCutsWithDefault()
      	exN03geom.SetMagField(magF)

      	print "Now geometry updated"

      	print position

#        gApplyUICommand("/vis/viewer/flush")
#        gApplyUICommand("/vis/scene/add/text 0 610 610 mm 20 0 0  " + "wxPython")

        gApplyUICommand("/gun/particle " +  str ( self.theParticle.GetStringSelection() ) )
        for i in self.processList:
#             print i, self.theProcesses.processCheck[i].GetValue()
             gProcessTable.SetProcessActivation(i, 1)
             if  self.theProcesses.processCheck[i].GetValue()  != True:
                 gProcessTable.SetProcessActivation(i, 0)
                 
        gApplyUICommand("/gun/energy  " + str ( self.energySpin.valAndUnit.GetValue() ) )

        eventNum = self.eventNo.slider.GetValue()
        for i in range(eventNum):
             pg.SetParticlePosition(G4ThreeVector(position, (i-eventNum/2)*5.*mm, 0.*cm))
             gRunManager.BeamOn(1)
#            sleep(0.01)
      	gApplyUICommand("/vis/viewer/update")

        
app = wx.PySimpleApp(False)
MyApp().Show()
app.MainLoop()
