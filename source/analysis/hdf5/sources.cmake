#------------------------------------------------------------------------------
# Module : G4hdf5
# Package: Geant4.src.G4analysis.G4hdf5
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4hdf5
  HEADERS
    G4Hdf5AnalysisManager.hh
    G4Hdf5AnalysisManager.icc
    G4Hdf5AnalysisReader.hh
    G4Hdf5AnalysisReader.icc
    G4Hdf5FileManager.hh
    G4Hdf5HnFileManager.hh
    G4Hdf5HnFileManager.icc
    G4Hdf5NtupleManager.hh
    G4Hdf5NtupleFileManager.hh
    G4Hdf5RNtupleManager.hh
    G4Hdf5RFileManager.hh
    g4hdf5_defs.hh
    g4hdf5.hh
  SOURCES
    G4Hdf5AnalysisManager.cc
    G4Hdf5AnalysisReader.cc
    G4Hdf5FileManager.cc
    G4Hdf5NtupleManager.cc
    G4Hdf5NtupleFileManager.cc
    G4Hdf5RFileManager.cc
    G4Hdf5RNtupleManager.cc
  GRANULAR_DEPENDENCIES
    G4globman
    G4intercoms
    G4analysismng
    G4hntools
  GLOBAL_DEPENDENCIES
    G4global
    G4intercoms
  LINK_LIBRARIES
    ${HDF5_LIBRARIES}
)

# List any source specific properties here
