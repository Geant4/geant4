#------------------------------------------------------------------------------
# sources.cmake
# Module : G4xml
# Package: Geant4.src.G4analysis.G4xml
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4xml
  HEADERS
    G4XmlAnalysisManager.hh
    G4XmlAnalysisManager.icc
    G4XmlAnalysisReader.hh
    G4XmlAnalysisReader.icc
    G4XmlFileManager.hh
    G4XmlHnFileManager.hh
    G4XmlHnFileManager.icc
    G4XmlNtupleFileManager.hh
    G4XmlNtupleManager.hh
    G4XmlRFileManager.hh
    G4XmlRNtupleManager.hh
    g4xml_defs.hh
    g4xml.hh
  SOURCES
    G4XmlAnalysisManager.cc
    G4XmlAnalysisReader.cc
    G4XmlNtupleFileManager.cc
    G4XmlFileManager.cc
    G4XmlNtupleManager.cc
    G4XmlRFileManager.cc
    G4XmlRNtupleManager.cc
  GRANULAR_DEPENDENCIES
    G4globman
    G4intercoms
    G4analysismng
    G4hntools
  GLOBAL_DEPENDENCIES
    G4global
    G4intercoms
  LINK_LIBRARIES
    ${EXPAT_LIBRARIES}
)

# List any source specific properties here
