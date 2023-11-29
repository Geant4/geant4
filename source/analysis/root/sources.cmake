# - G4root module build definition

# Define the Geant4 Module.
geant4_add_module(G4root
  PUBLIC_HEADERS
    G4RootAnalysisManager.hh
    G4RootAnalysisManager.icc
    G4RootAnalysisReader.hh
    G4RootAnalysisReader.icc
    G4RootFileDef.hh
    G4RootFileManager.hh
    G4RootHnFileManager.hh
    G4RootHnFileManager.icc
    G4RootHnRFileManager.hh
    G4RootHnRFileManager.icc
    G4RootMainNtupleManager.hh
    G4RootNtupleFileManager.hh
    G4RootNtupleManager.hh
    G4RootNtupleManager.icc
    G4RootPNtupleDescription.hh
    G4RootPNtupleManager.hh
    G4RootPNtupleManager.icc
    G4RootRFileDef.hh
    G4RootRFileManager.hh
    G4RootRNtupleManager.hh
    g4root_defs.hh
  SOURCES
    G4RootAnalysisManager.cc
    G4RootAnalysisReader.cc
    G4RootFileManager.cc
    G4RootMainNtupleManager.cc
    G4RootNtupleFileManager.cc
    G4RootNtupleManager.cc
    G4RootPNtupleManager.cc
    G4RootRFileManager.cc
    G4RootRNtupleManager.cc)

# Depends on ZLIB PUBLICally because of inclusion of tools/zlib in G4RootHnFileManager.icc
# (No, it't not obvious!)
geant4_module_link_libraries(G4root PUBLIC G4analysismng G4hntools G4globman G4tools ${ZLIB_LIBRARIES})
