# - G4visHepRep module build definition

# Define the Geant4 Module.
geant4_add_module(G4visHepRep
  PUBLIC_HEADERS
    G4HepRep.hh
    G4HepRepFile.hh
    G4HepRepFileSceneHandler.hh
    G4HepRepFileViewer.hh
    G4HepRepFileXMLWriter.hh
    G4HepRepMessenger.hh
    G4HepRepSceneHandler.hh
    G4HepRepViewer.hh
  SOURCES
    BHepRepWriter.cc
    DefaultHepRep.cc
    DefaultHepRepAction.cc
    DefaultHepRepAttDef.cc
    DefaultHepRepAttValue.cc
    DefaultHepRepAttribute.cc
    DefaultHepRepDefinition.cc
    DefaultHepRepFactory.cc
    DefaultHepRepInstance.cc
    DefaultHepRepInstanceTree.cc
    DefaultHepRepPoint.cc
    DefaultHepRepTreeID.cc
    DefaultHepRepType.cc
    DefaultHepRepTypeTree.cc
    DeflateOutputStreamBuffer.cc
    G4HepRep.cc
    G4HepRepFile.cc
    G4HepRepFileSceneHandler.cc
    G4HepRepFileViewer.cc
    G4HepRepFileXMLWriter.cc
    G4HepRepMessenger.cc
    G4HepRepSceneHandler.cc
    G4HepRepViewer.cc
    GZIPOutputStream.cc
    GZIPOutputStreamBuffer.cc
    IndentPrintWriter.cc
    XMLHepRepFactory.cc
    XMLHepRepWriter.cc
    XMLWriter.cc
    ZipOutputStream.cc
    ZipOutputStreamBuffer.cc)

geant4_module_link_libraries(G4visHepRep
  PUBLIC
    G4csg
    G4geometrymng
    G4materials
    G4modeling
    G4specsolids
    G4globman
    G4intercoms
    G4vis_management
    G4graphics_reps
  PRIVATE
    G4hepgeometry
    G4hits
    G4UIcommon
    G4tracking
    ${ZLIB_LIBRARIES})

# List any source specific properties here
