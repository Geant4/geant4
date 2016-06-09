#!/usr/bin/python

from Geant4 import *
import g4py.NISTmaterials
import g4py.ezgeom
from g4py.ezgeom import G4EzVolume
import g4py.ExN03pl
import g4py.ParticleGun
from time import *
import sys
from subprocess import *
import os

# ==================================================================
# intialize
# ==================================================================
def Configure():
  # ------------------------------------------------------------------
  # setup for materials
  # ------------------------------------------------------------------
  # simple materials for Qgeom
  g4py.NISTmaterials.Construct()

  # ------------------------------------------------------------------
  # setup for geometry
  # ------------------------------------------------------------------
  g4py.ezgeom.Construct()  # initialize

  # ------------------------------------------------------------------
  # setup for physics list of N03
  # ------------------------------------------------------------------
  g4py.ExN03pl.Construct()

  # ------------------------------------------------------------------
  # setup for primary generator action
  # ------------------------------------------------------------------
  g4py.ParticleGun.Construct()
  gControlExecute("gun.mac")

# ==================================================================
# constructing geometry
# ==================================================================
def ConstructGeom():

  print "* Constructing geometry..."
  # reset world material
  global absorber
  air= G4Material.GetMaterial("G4_AIR", 1)
  galactic = G4Material.GetMaterial("G4_Galactic", 1)
  absorber = {}  # material's dictionary to be used by a radiobutton
  aluminum = G4Material.GetMaterial("G4_Al", 1)
  iron   = G4Material.GetMaterial("G4_Fe", 1)
  silver = G4Material.GetMaterial("G4_Ag", 1)
  gold = G4Material.GetMaterial("G4_Au", 1)
  lead = G4Material.GetMaterial("G4_Pb", 1)
  water = G4Material.GetMaterial("G4_WATER", 1)
  absorber = {"air":air, "aluminum":aluminum, "iron":iron, "lead":lead, "water":water, "gold":gold
}
  g4py.ezgeom.SetWorldMaterial(galactic)
  g4py.ezgeom.ResizeWorld(120.*cm, 120.*cm, 100.*cm)
  # water phantom
  global water_phantom, water_phantom_pv

  water_phantom= G4EzVolume("WaterPhantom")
  water_phantom.CreateBoxVolume(water, 110.*cm, 110.*cm, 10.*cm)
  water_phantom_pv = water_phantom.PlaceIt(G4ThreeVector(0.,0.,0.*cm))

# ==================================================================
# main
# ==================================================================
# ------------------------------------------------------------------
# randum number
# ------------------------------------------------------------------
#print "Random numbers..."
rand_engine= Ranlux64Engine()
HepRandom.setTheEngine(rand_engine)
HepRandom.setTheSeed(20050830L)

# setup...

Configure()
ConstructGeom()

# ------------------------------------------------------------------
# go...
# ------------------------------------------------------------------
gRunManager.Initialize()



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
      pguide = pguide + "Parameter: " + str(iparam.GetParameterName())+ "   Type: " +           str(iparam.GetParameterType()) + '\n'
    guide = str(icommand.GetTitle()) + '\n' + pguide
    commandDic[command] = guide

  for i in range(1, ntree+1):
    itree= atree.GetTree(i)
    DumpTree(itree)

root_tree= gUImanager.GetTree()
DumpTree(root_tree)

# visualization
# OGLSX, VRML and HEPREP sceneHandlers are all created with names
gApplyUICommand("/vis/sceneHandler/create OGLSX OGLSX")
gApplyUICommand("/vis/sceneHandler/create VRML2FILE VRML")
gApplyUICommand("/vis/sceneHandler/create HepRepFile HEPREP")


#  OGLSX is the default so, viewer is created and volume is drawn
gApplyUICommand("/vis/viewer/create OGLSX oglsxviewer")
gApplyUICommand("/vis/viewer/select  oglsxviewer")
gApplyUICommand("/vis/drawVolume")
gApplyUICommand("/vis/scene/add/trajectories")

gApplyUICommand("/tracking/storeTrajectory 1")
gApplyUICommand("/vis/scene/endOfEventAction accumulate")
gApplyUICommand("/vis/scene/endOfRunAction accumulate")

gApplyUICommand("/vis/viewer/set/viewpointThetaPhi 90. 0.")
gApplyUICommand("/vis/viewer/set/style s")

# viewers VRML and Wired are tested by their envs vars
# if their envs var are set, then viewers are created and drawVolume

global heprepViewer, heprepDir, heprepName, vrmlViewer
heprepViewer = os.environ.get("G4HEPREPFILE_VIEWER")
heprepDir = os.environ.get("G4HEPREPFILE_DIR")
heprepName = os.environ.get("G4HEPREPFILE_NAME")
if heprepViewer is not None:
  gApplyUICommand("/vis/viewer/create HEPREP wired")
#  gApplyUICommand("/vis/drawVolume")

# VRML viewers name is user defined
vrmlDir = os.environ.get("G4VRML_DEST_DIR")
vrmlViewer = os.environ.get("G4VRMLFILE_VIEWER")
if vrmlViewer is not None:
  gApplyUICommand("/vis/viewer/create VRML vrmlviewer")
#  gApplyUICommand("/vis/drawVolume")

###############################
######  wxPython GUI ##########
###############################

import wx

class ComPanel(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent, -1)

    self.g4comText = wx.TextCtrl(parent, -1, "enter Command", pos=(10,10), size=(300, 30))

    self.g4comExec = wx.Button(parent, -1, "Execute", pos=(320,10), size=(60,30))
    self.g4comExec.Bind(wx.EVT_BUTTON, self.ExecuteCommand, self.g4comExec)

    self.comListBox = wx.ListBox(parent, -1, pos=(10,50), size=(300,200), choices=commandList, style=wx.LB_SINGLE)
    self.comListBox.SetSelection(1)
    self.comListBox.Bind(wx.EVT_LISTBOX, self.ShowGuide, self.comListBox)
    
    self.guide = wx.TextCtrl(parent, -1, "guidance", pos=(320,50), size=(300, 200), style =wx.TE_MULTILINE)
    self.guide.Bind(wx.EVT_LISTBOX, self.ShowGuide, self.guide)


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
    gApplyUICommand("/vis/viewer/pan " +  " -1. 0.  mm")
  def cmdLeft(self, event):
    gApplyUICommand("/vis/viewer/pan "   + " 1. 0. mm")

  def ViewerSelected(self, event):
    self.viewerName = event.GetString()

    if self.viewerName == "OpenGL":
      gApplyUICommand("/vis/viewer/select oglsxviewer")
      gApplyUICommand("/vis/drawVolume")
      gApplyUICommand("/vis/scene/add/trajectories")

      gApplyUICommand("/tracking/storeTrajectory 1")
      gApplyUICommand("/vis/scene/endOfEventAction accumulate")
      gApplyUICommand("/vis/scene/endOfRunAction accumulate")

    if self.viewerName == "VRML":
      gApplyUICommand("/vis/viewer/select vrmlviewer")

      gApplyUICommand("/vis/viewer/set/viewpointThetaPhi 90. 0.")
      gApplyUICommand("/vis/viewer/set/style s")
      gApplyUICommand("/vis/drawVolume")      
      gApplyUICommand("/vis/scene/add/trajectories")

      gApplyUICommand("/tracking/storeTrajectory 1")
      gApplyUICommand("/vis/scene/endOfEventAction accumulate")
      gApplyUICommand("/vis/scene/endOfRunAction accumulate")

    if self.viewerName == "Wired":

      gApplyUICommand("/vis/viewer/select wired")

      gApplyUICommand("/vis/viewer/set/viewpointThetaPhi 90. 0.")
      gApplyUICommand("/vis/viewer/set/style s")
      gApplyUICommand("/vis/drawVolume")
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
                             wx.DefaultSize, myList, 3, wx.RA_SPECIFY_ROWS)
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
#        self.theValue = str(self.intPart.GetValue()) + "." + str(self.manPart.GetValue()) + "  "
        self.theValue = float (self.intPart.GetValue()) + float(self.manPart.GetValue())/ 1000.
        self.theUnit =  self.unitSel.GetStringSelection()
        theText = "%.3f" % (self.theValue) + " " + self.theUnit
        self.valAndUnit.SetValue(theText)
# user may edit the Entry, so finnaly this value must be got
#        but it doesn't work for length but work for energy (gApplyUIcommnd)


# special class for this example to set/unset processes
class Processes(wx.Panel):
    def __init__(self, parent, myTitle, myList):
        wx.Panel.__init__(self, parent, -1)
        self.processCheck = {}
        self.sizer = wx.FlexGridSizer(rows=3)
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


############################
# WxPython Application class
############################

# main class to instantiate the above classes and pack them using nested sizers
        
class MyApp(wx.Frame):
    def __init__(self):
         wx.Frame.__init__(self, None, -1,  "Geant4Py", wx.DefaultPosition, size=(650,500))
         self.nb = wx.Notebook(self, -1, wx.DefaultPosition, size=(650,500),
                             style=
                             wx.NB_TOP # | wx.NB_MULTILINE
                             #wx.NB_BOTTOM
                             #wx.NB_LEFT
                             #wx.NB_RIGHT
                             )
         self.nb.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged)
         self.nb.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.OnPageChanging)
        
         panel = wx.Panel(self.nb)
	 self.nb.AddPage(panel, "Lesson1")

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

         self.runStart = wx.Button(panel, -1, "  Run Start", pos=(30,10), size=(80,30))
         self.Bind(wx.EVT_BUTTON, self.RunStart, self.runStart)
         bxsizer.Add(self.runStart, 0, wx.ALL)
# widgets
         

         materialList = absorber.keys()
         self.theMaterial = SelectOne(panel, "Materials", materialList)
  
         particleList = ['gamma', 'e-', 'proton']
         self.theParticle = SelectOne(panel, "Particles", particleList)

         self.processList = ["phot", "compt", "conv", "msc", "eIoni", "eBrem", "annihil", "hIoni"]
         self.theProcesses = Processes(panel, "Processes", self.processList)

         self.eventNo = Adjuster(panel, "Number of Events", 1 , 100 , 1)

         self.lengthUnit = {'micrometer':micrometer, 'mm':mm, 'cm':cm, 'm':m}
         self.thickSpin = FloatCounter(panel, "absorber thickness", self.lengthUnit.keys())

         self.energyUnit = { 'keV':keV, 'MeV':MeV, 'GeV':GeV, 'TeV':TeV, 'PeV':PeV}
         self.energySpin = FloatCounter(panel, "incident beam energy", self.energyUnit.keys())

# now sizers

         bxsizer.Add(wx.StaticLine(panel), 0, wx.EXPAND|wx.TOP|wx.BOTTOM, 5)
         bysizer.Add((10, -1))
         bysizer.Add(self.theMaterial, 0, wx.EXPAND, 10)
         bysizer.Add((10, -1))
         bysizer.Add(self.theParticle, 0, wx.EXPAND, 10)
         bysizer.Add((10, -1))         
         bysizer.Add(self.theProcesses.sizer, 0, wx.EXPAND, 10)
         
         bxsizer.Add(bysizer, 0, wx.EXPAND)

         bxsizer.Add(wx.StaticLine(panel), 0, wx.EXPAND|wx.TOP|wx.BOTTOM, 5)
         bxsizer.Add(self.energySpin.sizer, 0, wx.EXPAND)
         bxsizer.Add(wx.StaticLine(panel), 0, wx.EXPAND|wx.TOP|wx.BOTTOM, 5)         
         bxsizer.Add(self.thickSpin.sizer, 0, wx.EXPAND)         
         
         bxsizer.Add(wx.StaticLine(panel), 0, wx.EXPAND|wx.TOP|wx.BOTTOM, 5)
         bxsizer.Add(self.eventNo.sizer, 0, wx.EXPAND)

         self.solid = g4py.ezgeom.G4EzVolume.GetSold(water_phantom)
         
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


###### Run Action


    def RunStart(self, event):
         materialChosen =  str( self.theMaterial.GetStringSelection() ) 
         water_phantom.SetMaterial(absorber[materialChosen])
         if materialChosen == "water":
             water_phantom.SetColor(0., 0.9, 1.0)
         if materialChosen == "air":
             water_phantom.SetColor(0.9, 0.9, 1.0)
         if materialChosen == "lead":
             water_phantom.SetColor(0.2, 0.2, 0.2)
         if materialChosen == "iron":
             water_phantom.SetColor(0.7, 0.5, 0.7)
         if materialChosen == "aluminum":
             water_phantom.SetColor(.7, 0.9, 1.0)
         if materialChosen == "gold":
             water_phantom.SetColor(1., 0.9, .0)

         self.solid.SetZHalfLength(  self.thickSpin.theValue * self.lengthUnit[self.thickSpin.theUnit] / 2.0 )

#
         gApplyUICommand("/vis/scene/add/text 0 610 610 mm 20 0 0  " + "Geant4Py in Action")
         gApplyUICommand("/vis/viewer/flush")
         gApplyUICommand("/gun/particle " +  str ( self.theParticle.GetStringSelection() ) )
         for i in self.processList:
#             print i, self.theProcesses.processCheck[i].GetValue()
             gProcessTable.SetProcessActivation(i, 1)
             if  self.theProcesses.processCheck[i].GetValue()  != True:
                 gProcessTable.SetProcessActivation(i, 0)
                 
         gApplyUICommand("/gun/energy  " + str ( self.energySpin.valAndUnit.GetValue() ) )

         eventNum = self.eventNo.slider.GetValue()
         for i in range(eventNum):
             gunYZpos = str(i-eventNum/2) + ". -50. cm"
             gApplyUICommand("/gun/position 0. " + gunYZpos)
             gRunManager.BeamOn(1)
# next line is necessary for vrml to draw trajectories
         gApplyUICommand("/vis/viewer/update")

app = wx.PySimpleApp(False)
MyApp().Show()
app.MainLoop()
