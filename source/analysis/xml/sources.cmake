# - G4xml module build definition

# Define the Geant4 Module.
geant4_add_module(G4xml
  PUBLIC_HEADERS
    G4XmlAnalysisManager.hh
    G4XmlAnalysisManager.icc
    G4XmlAnalysisReader.hh
    G4XmlAnalysisReader.icc
    G4XmlFileManager.hh
    G4XmlHnFileManager.hh
    G4XmlHnFileManager.icc
    G4XmlHnRFileManager.hh
    G4XmlHnRFileManager.icc
    G4XmlNtupleFileManager.hh
    G4XmlNtupleManager.hh
    G4XmlRFileManager.hh
    G4XmlRFileManager.icc
    G4XmlRNtupleManager.hh
    g4xml_defs.hh
  SOURCES
    G4XmlAnalysisManager.cc
    G4XmlAnalysisReader.cc
    G4XmlNtupleFileManager.cc
    G4XmlFileManager.cc
    G4XmlNtupleManager.cc
    G4XmlRFileManager.cc
    G4XmlRNtupleManager.cc)

# EXPAT is a PUBLIC dependency by virtue of inclusion of tools/raxml
# (No, it's not obvious!) and exposure of that in G4Xml{RFileManager/AnalysisReader}
geant4_module_link_libraries(G4xml PUBLIC G4analysismng G4hntools G4globman G4tools ${EXPAT_LIBRARIES})
