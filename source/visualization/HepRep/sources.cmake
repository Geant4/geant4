#-----------------------------------------------------------------------
# Module : G4visHepRep
# Package: Geant4.src.G4visualization.G4visHepRep
#-----------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4visHepRep
  HEADERS
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
    ZipOutputStreamBuffer.cc
  GRANULAR_DEPENDENCIES
    G4csg
    G4geometrymng
    G4globman
    G4graphics_reps
    G4hits
    G4intercoms
    G4materials
    G4modeling
    G4specsolids
    G4tracking
    G4vis_management
  GLOBAL_DEPENDENCIES
    G4digits_hits
    G4geometry
    G4global
    G4graphics_reps
    G4intercoms
    G4materials
    G4modeling
    G4tracking
    G4vis_management
  LINK_LIBRARIES
    ${ZLIB_LIBRARIES}
  )

# List any source specific properties here
