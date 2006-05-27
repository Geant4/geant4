#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#   ExN03 : geant4/examples/novice/N03
#        using site-module packages
# ==================================================================
from Geant4 import *
import Qmaterials, NISTmaterials
import Qgeom,  ExN03geom
import ExN03pl, EMSTDpl
import ParticleGun, MedicalBeam
import sys
from time import *
from subprocess import *

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
#NISTmaterials.Construct()

# ------------------------------------------------------------------
# setup for geometry
# ------------------------------------------------------------------
# normal way for constructing user geometry


exN03geom= ExN03geom.ExN03DetectorConstruction()
gRunManager.SetUserInitialization(exN03geom)

# 2nd way, short-cut way

#ExN01geom.Construct()
#ExN03geom.Construct()

# magnetic field
#exN03geom.SetMagField(0.1 * tesla)

# ------------------------------------------------------------------
# setup for physics list
# ------------------------------------------------------------------
# normal way for constructing user physics list
exN03PL= ExN03pl.ExN03PhysicsList()
gRunManager.SetUserInitialization(exN03PL)

# 2nd way, short-cut way
#ExN01pl.Construct()
#EMSTDpl.Construct()

# ------------------------------------------------------------------
# setup for primary generator action
# ------------------------------------------------------------------
# normal way for constructing user physics list
#pgPGA= ParticleGun.ParticleGunAction()
#gRunManager.SetUserAction(pgPGA)
#pg= pgPGA.GetParticleGun()

# 2nd way, short-cut way
pg= ParticleGun.Construct()

# set parameters of particle gun
pg.SetParticleByName("e-")
pg.SetParticleEnergy(200.*MeV)
pg.SetParticlePosition(G4ThreeVector(-40.,0.,0.)*cm)
pg.SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.))

# medical beam
#beam= MedicalBeam.Construct()

# ------------------------------------------------------------------
# go...
# ------------------------------------------------------------------
gRunManager.Initialize()

# visualization
gApplyUICommand("/control/execute heprep.mac")

# beamOn
#gRunManager.BeamOn(3)


#TEST
#gProcessTable.SetProcessActivation("msc", 0)
#gProcessTable.SetProcessActivation("conv", 0)
#gProcessTable.SetProcessActivation("eBrem", 0)
#gProcessTable.SetProcessActivation("eIoni", 0)
#gProcessTable.SetProcessActivation("annihil", 0)


# creating widgets using grid layout

from Tkinter import *

class App(Frame):

  
  
  def init(self):

#title and header    row=0, 1
    title = Label(self, text="Geant4Py for Education \n@ H. Yoshida")
    title.grid(row=0, column=1, columnspan=3)
    header = Label(self, text="exampleN03 with Wired")
    header.grid(row=1, column=1, columnspan=3)
# number of layers

#absorber material selection row=2
    materialLabel = Label(self, bg="green", text="Absorber Material")
    materialLabel.grid(row=2, column=0, sticky=W)
    self.materialVar = StringVar()
    self.materialVar.set("Lead")
    ra1 = { }
    pos=1
    for i in ("Aluminium", "Lead"):
      ra1[i] = Radiobutton(self, text=i, variable=self.materialVar, value=i)
      ra1[i].grid(row=2, column=pos, sticky=W)
      pos=pos+1

#absorber thickness row=3
    thickLabel = Label(self, bg="green", text="Thickness (mm)")
    self.thickVar = DoubleVar()
    self.thickVar.set(1.0)
    thick = Scale(self, orient=HORIZONTAL, length=400, from_=0., to=100., resolution=0.05, tickinterval=10.0, digits=4, variable=self.thickVar)
    thickLabel.grid(row=3, column=0, sticky=W)
    thick.grid(row=3, column=1, columnspan=5, sticky=W)


#particle row=4
    particleLabel = Label(self, bg="green",  text="Particle")
    particleLabel.grid(row=4, column=0, sticky=W)
    self.particleVar = StringVar()
    self.particleVar.set("gamma")
    ra1 = { }
    pos1=1
    for i in ("gamma", "e-",  "e+", "mu-", "mu+"):
      ra1[i] = Radiobutton(self, text=i, variable=self.particleVar, value=i)
      ra1[i].grid(row=4, column=pos1, sticky=W)
      pos1=pos1+1

#energy row=5
    energyLabel = Label(self, bg="green",  text="Energy (MeV)")

    self.energyVar=StringVar()
    self.energyVar.set(1)
    energy = Scale(self, orient=HORIZONTAL, length=400, from_=0., to=1000., tickinterval=100.0, resolution=0.1, variable=self.energyVar, digits=5 )
    energyLabel.grid(row=5, column=0, sticky=W)
    energy.grid(row=5, column=1, columnspan=5, sticky=W)

#number of event row=6 
#    eventLabel = Label(self, bg="green",  text="One event")
#    self.eventVar=IntVar()
#    self.eventVar.set(3)
#    event = Scale(self,  orient=HORIZONTAL, length=400, from_=0, to=100, tickinterval=10, resolution=1, variable=self.eventVar )
#    eventLabel.grid(row=6, column=0, sticky=W)
#    event.grid(row=6, column=1, columnspan=5, sticky=W)

#start a run button row=0
    startBut = Button(self, bg="orange", text="Start a run", command=self.cmd_beamOn)
    startBut.grid(row=0, column=0, sticky=W)

#Zoom in/out Pan X Y row=13
    visLabel = Label(self, text="Wired", bg="orange")
    launchBut = Button(self, text="Launch", bg="orange", command=self.cmd_wired)
#    shrinkBut = Button(self, text="Zoom out", command=self.cmd_shrink)
    visLabel.grid(row=13, column=0, sticky=W)
    launchBut.grid(row=13, column=1, sticky=W)
#    shrinkBut.grid(row=13, column=2, sticky=W)
#    panLabel = Label(self, text="Pan X Y(mm)")
#    self.panXYVar = StringVar()
#    panXYEnt = Entry(self, textvariable=self.panXYVar, width=6)
#    panBut = Button(self, bg="orange", text="OK", command=self.cmd_pan)
#    panLabel.grid(row=13, column=3, sticky=W)
#    panXYEnt.grid(row=13, column=4)
#    panBut.grid(row=13, column=5)
#Geant4 command entry row = 14
    g4comLabel = Label(self, text="Geant4 command", bg="orange")
    self.g4commandVar = StringVar()
    commandEntry = Entry(self, textvariable=self.g4commandVar, width=15)
    comBut = Button(self, bg="orange", text="Execute", command=self.cmd_g4command)
    g4comLabel.grid(row=14, column=0, sticky=W)
    commandEntry.grid(row=14, column=1, columnspan=3, sticky=E+W)
    comBut.grid(row=14, column=5)

# process activate row 8
    processLabel=Label(self, text="Process on/off", bg="green")
    processLabel.grid(row=8, column=0, sticky=W)
    procTab = {}
    
    self.processList = ["phot", "compt", "conv", "msc", "eIoni", "eBrem", "annihil","muIoni", "muBrems"]
    pos=1
    self.processVar = {}
    for i in self.processList:
      self.processVar[i] = IntVar()
      procTab[i] = Checkbutton(self, text=i, variable=self.processVar[i], command=self.cmd_setProcess)
      if pos <= 3:
        procTab[i].grid(row=8, column=pos, sticky=W)
      if 4<= pos <= 7:
        procTab[i].grid(row=9, column=pos-3, sticky=W)
      if pos >= 8:
        procTab[i].grid(row=10, column=pos-7, sticky=W)
      pos=pos+1
      procTab[i].select()
# set cuts row 11
    cutLabel = Label(self, bg="green",  text="Cut (mm)")

    self.cutVar=DoubleVar()
    self.cutVar.set(1.)
    cut = Scale(self, orient=HORIZONTAL, length=400, from_=0., to=10., tickinterval=1., resolution=0.005, variable=self.cutVar, digits=5 )
    cutLabel.grid(row=11, column=0, sticky=W)
    cut.grid(row=11, column=1, columnspan=5, sticky=W)

# set mag field row 12
    magLabel = Label(self, bg="green",  text="Magnetic (T)")

    self.magVar=DoubleVar()
    self.magVar.set(0.)
    mag = Scale(self, orient=HORIZONTAL, length=400, from_=0., to=5., tickinterval=1., resolution=0.1, variable=self.magVar, digits=3 )
    magLabel.grid(row=12, column=0, sticky=W)
    mag.grid(row=12, column=1, columnspan=5, sticky=W)


#exit row = 0    
    exitBut = Button(self, bg="red", text="End all", command=sys.exit)
    exitBut.grid(row=0, column=5, sticky=W)

#on Run butto do...    
  def cmd_beamOn(self):
      exN03geom.SetAbsorberMaterial(self.materialVar.get())
      exN03geom.SetAbsorberThickness(self.thickVar.get()  * mm/2.0)
      exN03geom.UpdateGeometry()
      exN03PL.SetDefaultCutValue(self.cutVar.get() * mm)
      exN03PL.SetCutsWithDefault()
      exN03geom.SetMagField(self.magVar.get() * tesla)
      print "Now geometry updated"

#      gApplyUICommand("/vis/viewer/flush")
      self.cmd_particle(self.particleVar.get())
      self.cmd_energy(self.energyVar.get())
#      gApplyUICommand("/vis/scene/add/text 0 610 610 mm 20 0 0  " + self.materialVar.get() + " = " + str(self.thickVar.get()) + "mm " + self.particleVar.get() + " = "+self.energyVar.get() + "MeV")

#      eventNum = self.eventVar.get()
#      for i in range(eventNum):

#        gunYZpos = str((i-eventNum/2)*5) + ". 0. mm"
#        gApplyUICommand("/gun/position -40. " + gunYZpos)
      gRunManager.BeamOn(1)
      gApplyUICommand("/vis/viewer/flush") # test HepRepFile
#        sleep(0.01)

  def cmd_setProcess(self):
    for i in self.processList:
      if self.processVar[i].get() == 0:
         gProcessTable.SetProcessActivation(i, 0)
#         print "Process " + i + " inactivated"
      else:
         gProcessTable.SetProcessActivation(i, 1)
#         print "Process " + i + " activated"
      gApplyUICommand("/run/physicsModified")
      print "Physics modified"

      
  def cmd_g4command(self):
    gApplyUICommand(self.g4commandVar.get())
      
  def cmd_particle(self, particle):
    gApplyUICommand("/gun/particle " + particle)


  def cmd_energy(self, penergy):
    gApplyUICommand("/gun/energy " + penergy + " MeV")


  def cmd_wired(self):
    Popen("/home/yoshidah/Wired/bin/wired" +  " -file g4hr_00.wir.heprep", shell=True) 

#  def cmd_pan(self):
#    gApplyUICommand("/vis/viewer/pan " + self.panXYVar.get() + " "  + " mm")


#  def cmd_shrink(self):
#    gApplyUICommand("/vis/viewer/zoom 0.8")


    
  def __init__(self, master=None):
    Frame.__init__(self, master)
    self.init()
    self.grid()

    
app = App()
app.mainloop()
