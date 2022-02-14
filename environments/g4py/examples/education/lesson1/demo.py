#!/usr/bin/python
# ==================================================================
# python script for "measurement" of mass attenuation coefficient
#
#   
#   -  using site-module packages
# ==================================================================
from Geant4 import *
import g4py.NISTmaterials
import g4py.ezgeom
from g4py.ezgeom import G4EzVolume
import g4py.EMSTDpl
import g4py.ParticleGun
from time import *
import sys

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
  #g4py.Qgeom.Construct()
  g4py.ezgeom.Construct()  # initialize

  # ------------------------------------------------------------------
  # setup for physics list
  # ------------------------------------------------------------------
  g4py.EMSTDpl.Construct()

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
  absorber = {"air":air, "aluminum":aluminum, "iron":iron, "lead":lead, "water":water, "gold":gold}
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
print "Random numbers..."
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

# visualization not here but after "Start a run" button 
gControlExecute("oglx.mac")
#gControlExecute("vrml.mac")

# creating widgets using grid layout

from Tkinter import *

class App(Frame):
  
  
  def init(self):

#title and header    row=0, 1
    title = Label(self, text="Geant4Py for Education @ H. Yoshida Naruto Univ. of Education")
    title.grid(row=0, column=1, columnspan=5)
    header = Label(self, text="Measurement of Mass Attenuation Coefficient")
    header.grid(row=1, column=1, columnspan=5)

#material selection row=2
    materialLabel = Label(self, bg="green", text="Material")
    materialLabel.grid(row=2, column=0, sticky=W)
    self.materialVar = StringVar()
    self.materialVar.set("water")
    ra1 = { }
    pos=1
    for i in absorber.keys():
      ra1[i] = Radiobutton(self, text=i, variable=self.materialVar, value=i)
      ra1[i].grid(row=2, column=pos, sticky=W)
      pos=pos+1

#absorber thickness row=3
    thickLabel = Label(self, bg="green", text="Thickness (mm)")
    self.thickVar = DoubleVar()
    self.thickVar.set(100.0)
    thick = Scale(self, orient=HORIZONTAL, length=400, from_=0., to=100., resolution=0.05, tickinterval=10.0, digits=4, variable=self.thickVar)
    thickLabel.grid(row=3, column=0, sticky=W)
    thick.grid(row=3, column=1, columnspan=5, sticky=W)
    
#get logical volume and set its half length
    self.solid = g4py.ezgeom.G4EzVolume.GetSold(water_phantom)

#particle row=4
    particleLabel = Label(self, bg="green",  text="Particle")
    particleLabel.grid(row=4, column=0, sticky=W)
    self.particleVar = StringVar()
    self.particleVar.set("gamma")
    ra1 = { }
    pos1=1
    for i in ("gamma", "e-"):
      ra1[i] = Radiobutton(self, text=i, variable=self.particleVar, value=i)
      ra1[i].grid(row=4, column=pos1, sticky=W)
      pos1=pos1+1

#energy row=5
    energyLabel = Label(self, bg="green",  text="Energy (MeV)")

    self.energyVar=StringVar()
    self.energyVar.set(1)
    energy = Scale(self, orient=HORIZONTAL, length=400, from_=0., to=100., tickinterval=10.0, resolution=0.1, variable=self.energyVar, digits=4 )
    energyLabel.grid(row=5, column=0, sticky=W)
    energy.grid(row=5, column=1, columnspan=5, sticky=W)

#number of event row=6 
    eventLabel = Label(self, bg="green",  text="Events")
    self.eventVar=IntVar()
    event = Scale(self,  orient=HORIZONTAL, length=400, from_=1, to=100, tickinterval=10, resolution=1, variable=self.eventVar )
    eventLabel.grid(row=6, column=0, sticky=W)
    event.grid(row=6, column=1, columnspan=5, sticky=W)

#start a run button row=7
    startBut = Button(self, bg="orange", text="Start a run", command=self.cmd_beamOn)
    startBut.grid(row=0, column=0, sticky=W)

#Zoom in/out Pan X Y row=8
    visLabel = Label(self, text="viewer", bg="orange")
    expandBut = Button(self, text="Zoom in", command=self.cmd_expand)
    shrinkBut = Button(self, text="Zoom out", command=self.cmd_shrink)
    visLabel.grid(row=8, column=0, sticky=W)
    expandBut.grid(row=8, column=1, sticky=W)
    shrinkBut.grid(row=8, column=2, sticky=W)

    upBut = Button(self, text="Up", command=self.cmd_up)
    downBut = Button(self, text="Down", command=self.cmd_down)
    upBut.grid(row=8, column=3, sticky=W)
    downBut.grid(row=8, column=4, sticky=W)

    leftBut = Button(self, text="Left", command=self.cmd_left)
    rightBut = Button(self, text="Right", command=self.cmd_right)
    leftBut.grid(row=8, column=5, sticky=W)
    rightBut.grid(row=8, column=6, sticky=W)
# later
#    resetBut = Button(self, text="Reset", command=self.cmd_reset)
#    resetBut.grid(row=8, column=7, sticky=W)


#    panLabel = Label(self, text="Pan X Y (mm)")
#    self.panXYVar = StringVar()
#    panXYEnt = Entry(self, textvariable=self.panXYVar)
#    panBut = Button(self, bg="orange", text="OK", command=self.cmd_pan)
#    panLabel.grid(row=8, column=3, sticky=W)
#    panXYEnt.grid(row=8, column=4)
#    panBut.grid(row=8, column=5)
#Geant4 command entry row = 9
#    g4comLabel = Label(self, text="Geant4 command")
#    self.g4commandVar = StringVar()
#    commandEntry = Entry(self, textvariable=self.g4commandVar)
#    comBut = Button(self, bg="orange", text="Execute", command=self.cmd_g4command)
#    g4comLabel.grid(row=9, column=0, sticky=W)
#    commandEntry.grid(row=9, column=1, columnspan=4, sticky=E+W)
#    comBut.grid(row=9, column=5)

#exit row = 10    
    exitBut = Button(self, bg="red", text="End all", command=sys.exit)
    exitBut.grid(row=0, column=6, sticky=W)

#on Run butto do...    
  def cmd_beamOn(self):
      materialChosen = self.materialVar.get()
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
          
      self.solid.SetZHalfLength(self.thickVar.get() * mm/2.0)
#      gControlExecute("oglx.mac") #draw for each run
      gApplyUICommand("/vis/viewer/flush")

      self.cmd_particle(self.particleVar.get())
      self.cmd_energy(self.energyVar.get())
# TODO later to reflesh text
      gApplyUICommand("/vis/scene/add/text 0 610 610 mm 20 0 0  " + "                                            ")
      gApplyUICommand("/vis/scene/add/text 0 610 610 mm 20 0 0  " + self.materialVar.get() + " = " + str(self.thickVar.get()) + "mm " + self.particleVar.get() + " = "+self.energyVar.get() + "MeV")

      eventNum = self.eventVar.get()
      for i in range(eventNum):
        gunYZpos = str(i-eventNum/2) + ". -20. cm"
        gApplyUICommand("/gun/position 0. " + gunYZpos)
        gRunManager.BeamOn(1)
        sleep(0.01)
#      self.cmd_expand() #Zoom in to the last diaplayed OGLSX
#      self.cmd_shrink()


      
  def cmd_g4command(self):
    gApplyUICommand(self.g4commandVar.get())
      
  def cmd_particle(self, particle):
    gApplyUICommand("/gun/particle " + particle)


  def cmd_energy(self, penergy):
    gApplyUICommand("/gun/energy " + penergy + " MeV")


  def cmd_expand(self):
    gApplyUICommand("/vis/viewer/zoom 1.2")

  def cmd_up(self):
    gApplyUICommand("/vis/viewer/pan "   + " 0.  10. mm")

  def cmd_down(self):
    gApplyUICommand("/vis/viewer/pan " +  " 0. -10.  mm")

  def cmd_right(self):
    gApplyUICommand("/vis/viewer/pan " +  " -1. 0.  mm")

  def cmd_left(self):
    gApplyUICommand("/vis/viewer/pan "   + " 1. 0. mm")


  def cmd_shrink(self):
    gApplyUICommand("/vis/viewer/zoom 0.8")


#  def cmd_reset(self):
#    gApplyUICommand("/vis/viewer/pan "   + " 0. 0. mm")

    
  def __init__(self, master=None):
    Frame.__init__(self, master)
    self.init()
    self.grid()

    
app = App()
app.mainloop()
