# $Id: __init__.py,v 1.1 2006-02-27 09:56:05 kmura Exp $
"""
# ==================================================================
#  [Geant4] module package
#
#  Welcome to Geant4Py.
#
#  This package contains a set of Python interface with Geant4.
#  Hava A Fun!
# ==================================================================
"""

# ==================================================================
# docs
# ==================================================================
__version__ = '0.3.1'
__date__    = '19/Jan/2006'
__author__  = 'K.Murakami (Koichi.Murakami@kek.jp)'

print """
=============================================================
  Welcome to Geant4Py (Geant4 Python Interface)
  
  Version : %s
  Date    : %s
  Author  : %s
=============================================================
""" % ( __version__, __date__, __author__ )

# ==================================================================
# import submodules
# ==================================================================
from G4interface      import *
from G4global         import *
from G4run            import *
from G4event          import *
from G4tracking       import *
from G4track          import *
from G4particles      import *
from G4processes      import *
from G4geometry       import *
from G4materials      import *
from G4digits_hits    import *
from G4visualization  import *
from G4graphics_reps  import *
from HEPUnit          import *
from colortable       import *
from g4viscp          import *

# ==================================================================
# globals, which start with "g"
# ==================================================================
# ------------------------------------------------------------------
# gRunManager
# ------------------------------------------------------------------
global gRunManager
if(G4RunManager.GetRunManager() == None):
  gRunManager= G4RunManager()
else:
  gRunManager= G4RunManager.GetRunManager()

# ------------------------------------------------------------------
# gEventManager
# ------------------------------------------------------------------
global gEventManager
gEventManager= G4EventManager.GetEventManager()

# ------------------------------------------------------------------
# gStackManager
# ------------------------------------------------------------------
global gStackManager
gStackManager= gEventManager.GetStackManager()

# ------------------------------------------------------------------
# gTrackingManager
# ------------------------------------------------------------------
global gTrackingManager
gTrackingManager= gEventManager.GetTrackingManager()

# ------------------------------------------------------------------
# gStateManager
# ------------------------------------------------------------------
global gStateManager
gStateManager= G4StateManager.GetStateManager()

# ------------------------------------------------------------------
# gTransportationManager
# ------------------------------------------------------------------
global gTransportationManager
gTransportationManager= G4TransportationManager.GetTransportationManager()

# ------------------------------------------------------------------
# gParticleTable
# ------------------------------------------------------------------
global gParticleTable
gParticleTable= G4ParticleTable.GetParticleTable()

# ------------------------------------------------------------------
# gProcessTable
# ------------------------------------------------------------------
global gProcessTable
gProcessTable= G4ProcessTable.GetProcessTable()

# ------------------------------------------------------------------
# gNistManager (since 7.1)
# ------------------------------------------------------------------
global gNistManager

material_class_list= dir(G4materials)
qfind= (material_class_list.count("G4NistManager") >0)
if(qfind) :
  gNistManager= G4NistManager.Instance()


# ------------------------------------------------------------------
# gVisManager
# ------------------------------------------------------------------
global gVisManager
global opengl_ix, opengl_sx, opengl_ixm, opengl_sxm, vrml1, vrml2, dawn, \
       heprep_xml, heprep_file, atree, raytracer, raytracer_x

visdriver_list = dir(G4visualization)
q_opengl_ix   = (visdriver_list.count("G4OpenGLImmediateX") > 0)
q_opengl_sx   = (visdriver_list.count("G4OpenGLStoredX") >0)
q_opengl_ixm  = (visdriver_list.count("G4OpenGLImmediateXm") > 0)
q_opengl_sxm  = (visdriver_list.count("G4OpenGLStoredXm") >0)
q_raytracer_x = (visdriver_list.count("G4RayTracerX") >0)
  
if(G4VisManager.GetConcreteInstance() == None):
  gVisManager= G4VisManager()
  if(q_opengl_ix):
    opengl_ix= G4OpenGLImmediateX()
  if(q_opengl_sx):
    opengl_sx= G4OpenGLStoredX()
  if(q_opengl_ixm):
    opengl_ixm= G4OpenGLImmediateXm()
  if(q_opengl_sxm):
    opengl_sxm= G4OpenGLStoredXm()
  if(q_raytracer_x):
    raytracer_x= G4RayTracerX()

  vrml1= G4VRML1File()
  vrml2= G4VRML2File()
  dawn= G4DAWNFILE()
  heprep_xml= G4HepRep()
  heprep_file= G4HepRepFile()
  atree= G4ASCIITree()
  raytracer= G4RayTracer()

  if(q_opengl_ix):
    gVisManager.RegisterGraphicsSystem(opengl_ix)
  if(q_opengl_sx):
    gVisManager.RegisterGraphicsSystem(opengl_sx)
  if(q_opengl_ixm):
    gVisManager.RegisterGraphicsSystem(opengl_ixm)
  if(q_opengl_sxm):
    gVisManager.RegisterGraphicsSystem(opengl_sxm)
  if(q_raytracer_x):
    gVisManager.RegisterGraphicsSystem(raytracer_x)

  gVisManager.RegisterGraphicsSystem(vrml1)
  gVisManager.RegisterGraphicsSystem(vrml2)
  gVisManager.RegisterGraphicsSystem(dawn)
  gVisManager.RegisterGraphicsSystem(heprep_xml)
  gVisManager.RegisterGraphicsSystem(heprep_file)
  gVisManager.RegisterGraphicsSystem(atree)
  gVisManager.RegisterGraphicsSystem(raytracer)
    
  gVisManager.Initialize()

# ------------------------------------------------------------------
# functions
# ------------------------------------------------------------------
global gApplyUICommand
gApplyUICommand= G4interface.ApplyUICommand

global gGetCurrentValues
gGetCurrentValues= G4interface.GetCurrentValues

global gStartUISession
gStartUISession= G4interface.StartUISession

# ==================================================================
# extentions
# ==================================================================

# ------------------------------------------------------------------
# generate one event
# ------------------------------------------------------------------
def OneEvent(self):
  "generate one event."
  self.BeamOn(1)

G4RunManager.OneEvent = OneEvent

# ------------------------------------------------------------------
# list material information
# ------------------------------------------------------------------
global gMaterialTable
gMaterialTable= G4Material.GetMaterialTable()

global gElementTable
gElementTable= G4Element.GetElementTable()

def ListMaterial(self):
  "list materials."
  n_materials= len(gMaterialTable)
  print " +------------------------------------------------------------------"
  print " |       Table of G4Material-s (%d materails defined)" % (n_materials)
  for i in range(0, n_materials) :
    material= gMaterialTable[i]
    print " |--------------------------------------------------------"\
          "----------"
    print " | %s: %s" % (material.GetName(),
                         G4BestUnit(material.GetDensity(),"Volumic Mass"))

    elementVec= material.GetElementVector()
    fractionVec= material.GetFractionVector()
    abundanceVec= material.GetVecNbOfAtomsPerVolume()
    totNAtoms= material.GetTotNbOfAtomsPerVolume()

    n_elements= len(elementVec)
    for j in range(0, n_elements):
      print " | + (%1d) %s(%s): A=%4.1f, N=%5.1f, " \
            "Frac.=(%4.1f%%m,%4.1f%%a)" % \
            (j+1, elementVec[j].GetName(), elementVec[j].GetSymbol(),
             elementVec[j].GetZ(),
             elementVec[j].GetN(),
             fractionVec[j]/HEPUnit.perCent,
             abundanceVec[j]/totNAtoms/HEPUnit.perCent)

  print " +------------------------------------------------------------------"

G4MaterialTable.ListMaterial = ListMaterial

# ------------------------------------------------------------------
# execute G4 macro
# ------------------------------------------------------------------
def gControlExecute(g4mac):
  ui_command= "/control/execute " + g4mac
  gApplyUICommand(ui_command)
  
