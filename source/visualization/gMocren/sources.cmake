# - G4GMocren module definition

# Define the Geant4 Module.
geant4_add_module(G4GMocren
  PUBLIC_HEADERS
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
    G4GMocrenMessenger.cc)

geant4_module_link_libraries(G4GMocren
  PUBLIC
    G4hits
    G4intercoms
    G4vis_management
    G4geometrymng
    G4globman
  PRIVATE
    G4graphics_reps
    G4csg
    G4modeling
    G4materials
    G4navigation
    G4hepgeometry
    G4detutils
    G4tracking)
