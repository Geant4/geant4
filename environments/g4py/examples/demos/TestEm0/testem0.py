#!/usr/bin/python
# ==================================================================
# python script for TestEm0 python version
#
# ==================================================================
import Geant4 as g4
import TestEm0

# ==================================================================
# user actions in python
# ==================================================================


# ==================================================================
# main
# ==================================================================

myDC= TestEm0.DetectorConstruction()
g4.gRunManager.SetUserInitialization(myDC)

myPL= TestEm0.PhysicsList()
g4.gRunManager.SetUserInitialization(myPL)

# set user actions...
myPGA= TestEm0.PrimaryGeneratorAction(myDC)
g4.gRunManager.SetUserAction(myPGA)
        
myRA= TestEm0.RunAction(myDC,myPGA)

# set user action classes
g4.gRunManager.SetUserAction(myRA)

  
  
g4.gRunManager.Initialize()

pg = g4.G4ParticleGun()
    
materialList = TestEm0.getMaterialTable();

particleList = TestEm0.getParticleTable()

enrgyList = ["eV","keV","MeV","GeV","TeV","PeV"]

cutsList = ["um", "mm" ,  "cm",  "m",  "km"]




# GUI

from Tkinter import *
class App(Frame):

  
  def init(self):

# title and header
    title = Label(self, text="TestEm0 empowered by Geant4Py\n\n\n")
    title.grid(row=0, column=1, columnspan = 4)

# particle list box 
    particle_title =  Label(self, text="Particle")
    particle_title.grid(row=2, column=0)

    particleFrame = Frame(self)
    scrollbar2 = Scrollbar(particleFrame)
    scrollbar2.pack(side = RIGHT,  fill = Y)
    self.particleListBox = Listbox(particleFrame, yscrollcommand=scrollbar2.set, exportselection=FALSE,height = 6)
    self.particleListBox.pack(side = LEFT)
    for item in particleList:
       self.particleListBox.insert(END, item)
    scrollbar2.config(command=self.particleListBox.yview)
    particleFrame.grid(row=3, column=0)
    self.particleListBox.select_set(0)

# separator frame
    fblank = Frame(self,width = 40)
    fblank.grid(row=3,column=1)

# material list box 
    detmaterial_title =  Label(self, text="Material")
    detmaterial_title.grid(row=2, column=2)

    materialFrame = Frame(self)
    scrollbar = Scrollbar(materialFrame)
    scrollbar.pack(side = RIGHT, fill = Y)
    self.materialListBox = Listbox(materialFrame,  yscrollcommand=scrollbar.set, exportselection=FALSE, height = 6)
    self.materialListBox.pack(side = LEFT, fill = Y)
    for item in materialList:
       self.materialListBox.insert(END, item)
    scrollbar.config(command=self.materialListBox.yview)
    materialFrame.grid(row=3, column=2)
    self.materialListBox.select_set(0)

# separator frame
    fblank = Frame(self,width = 40)
    fblank.grid(row=3,column=3)

# energy
    fEnergy = Frame(self)
    energyLabel = Label(self,  text="Energy")
    energyLabel.grid(row = 2, column = 4)

    scrollbarEnergy = Scrollbar(fEnergy)
    scrollbarEnergy.pack(side = RIGHT,  fill = Y)
    self.energyEntry = Entry(fEnergy, width=  8 );
    self.energyEntry.pack(side = TOP)
    self.energyEntry.insert(0, "1.0")

    self.energyListBox = Listbox(fEnergy,  yscrollcommand=scrollbarEnergy.set,exportselection=FALSE,width=8,height = 5)
    self.energyListBox.pack(side = BOTTOM  )
    for item in enrgyList:
        self.energyListBox.insert(END, item)
    scrollbarEnergy.config(command=self.energyListBox.yview)
    fEnergy.grid(row = 3, column = 4 )
    self.energyListBox.select_set(0)

# separator frame
    fblank = Frame(self,width = 40)
    fblank.grid(row=3,column=5)

# cuts 
    fCuts = Frame(self)
    cutsLabel = Label(self,  text="Cuts",  width=  8)
    cutsLabel.grid(row = 2, column = 6)

    scrollbarCuts = Scrollbar(fCuts)
    scrollbarCuts.pack(side = RIGHT,  fill = Y)
    self.cutsEntry = Entry(fCuts, width=  8);
    self.cutsEntry.pack(side = TOP)
    self.cutsEntry.insert(0, "1.0")

    self.cutsListBox = Listbox(fCuts,  width=  8 ,yscrollcommand=scrollbarCuts.set,exportselection=FALSE,height = 5)
    self.cutsListBox.pack(side = BOTTOM  )
    for item in cutsList:
        self.cutsListBox.insert(END, item)
    scrollbarCuts.config(command=self.cutsListBox.yview)
    fCuts.grid(row = 3, column = 6 )
    self.cutsListBox.select_set(0)

# separator frame
    fblank = Frame(self,height = 40)
    fblank.grid(row=4,column=0)

# start a run button
    startBut = Button(self, bg="green", text="Start a run", command=self.cmd_beamOn)
    startBut.grid(row=5, column=2, sticky=W)
    
# exit button
    exitBut = Button(self, bg="grey", text="Exit", command=self.quit)
    exitBut.grid(row=5,  column=6,  sticky=E)
    
  def __init__(self, master=None):
    Frame.__init__(self, master)
    self.init()
    self.grid()

  def cmd_beamOn(self):
    
    # get and set particle
    if self.particleListBox.curselection():
        index =int(self.particleListBox.curselection()[0])
        g4.gApplyUICommand("/gun/particle  " +  particleList[index])

    # get and set detector Material
    if self.materialListBox.curselection():
        index =int(self.materialListBox.curselection()[0])
        g4.gApplyUICommand("/testem/det/setMat  " +  materialList[index])
	
    # get and set energy
    energy = self.energyEntry.get()
    if self.energyListBox.curselection():
        index = int(self.energyListBox.curselection()[0])
        unity =   enrgyList[index]
        g4.gApplyUICommand("/gun/energy " + energy + " " + unity)

    # get and set cuts
    cuts = self.cutsEntry.get()
    if self.cutsListBox.curselection():
        index = int(self.cutsListBox.curselection()[0])
        unity =   cutsList[index]
        g4.gApplyUICommand("/testem/phys/setCuts " + cuts + " " + unity)

    # run beamOn
    g4.gRunManager.BeamOn(1)


app = App()
app.mainloop()

