#------------------------------------------------------------------------------
# Module : G4GMocren
# Package: Geant4.src.G4visualization.G4GMocren
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4GMocren
  HEADERS
    G4GMocrenFile.hh
    G4GMocrenFileCTtoDensityMap.hh
    G4GMocrenFileSceneHandler.hh
    G4GMocrenFileViewer.hh
    G4GMocrenIO.hh
    G4GMocrenMessenger.hh
    G4GMocrenTouchable.hh
  SOURCES
    G4GMocrenFile.cc
    G4GMocrenFileCTtoDensityMap.cc
    G4GMocrenFileSceneHandler.cc
    G4GMocrenFileViewer.cc
    G4GMocrenIO.cc
    G4GMocrenMessenger.cc
  GRANULAR_DEPENDENCIES
    G4csg
    G4detutils
    G4digits
    G4event
    G4geometrymng
    G4globman
    G4graphics_reps
    G4hepnumerics
    G4hits
    G4intercoms
    G4materials
    G4modeling
    G4navigation
    G4partman
    G4specsolids
    G4tracking
    G4vis_management
  GLOBAL_DEPENDENCIES
    G4digits_hits
    G4event
    G4geometry
    G4global
    G4graphics_reps
    G4intercoms
    G4materials
    G4modeling
    G4particles
    G4tracking
    G4vis_management
)

# List any source specific properties here

