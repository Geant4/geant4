#------------------------------------------------------------------------------
# Module : G4root
# Package: Geant4.src.G4analysis.G4root
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4root
  HEADERS
    G4RootAnalysisManager.hh
    G4RootAnalysisManager.icc
    G4RootAnalysisReader.hh
    G4RootAnalysisReader.icc
    G4RootFileManager.hh
    G4RootMainNtupleManager.hh
    G4RootNtupleManager.hh
    G4RootPNtupleDescription.hh
    G4RootPNtupleManager.hh
    G4RootRFileManager.hh
    G4RootRNtupleManager.hh
    g4root_defs.hh
    g4root.hh
  SOURCES
    G4RootAnalysisManager.cc
    G4RootAnalysisReader.cc
    G4RootFileManager.cc
    G4RootMainNtupleManager.cc
    G4RootNtupleManager.cc
    G4RootPNtupleManager.cc
    G4RootRFileManager.cc
    G4RootRNtupleManager.cc
  GRANULAR_DEPENDENCIES
    G4globman
    G4intercoms
    G4analysismng
    G4hntools
  GLOBAL_DEPENDENCIES
    G4global
    G4intercoms
  LINK_LIBRARIES
    ${ZLIB_LIBRARIES}
)

# List any source specific properties here
