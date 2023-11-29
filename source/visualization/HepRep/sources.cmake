# - G4visHepRep module build definition

# Define the Geant4 Module.
geant4_add_module(G4visHepRep
  PUBLIC_HEADERS
    G4HepRepFile.hh
    G4HepRepFileXMLWriter.hh
  PRIVATE_HEADERS
    G4HepRepFileSceneHandler.hh
    G4HepRepFileViewer.hh
    G4HepRepMessenger.hh
  SOURCES
    G4HepRepFile.cc
    G4HepRepFileSceneHandler.cc
    G4HepRepFileViewer.cc
    G4HepRepFileXMLWriter.cc
    G4HepRepMessenger.cc)

geant4_module_link_libraries(G4visHepRep
  PUBLIC
    G4globman
    G4vis_management
  PRIVATE
    G4csg
    G4geometrymng
    G4graphics_reps
    G4hits
    G4intercoms
    G4materials
    G4modeling
    G4specsolids
    G4tracking)

# List any source specific properties here
