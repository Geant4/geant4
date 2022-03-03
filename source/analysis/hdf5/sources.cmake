# - G4hdf5 module build definition

# Define the Geant4 Module.
geant4_add_module(G4hdf5
  PUBLIC_HEADERS
    G4Hdf5AnalysisManager.hh
    G4Hdf5AnalysisManager.icc
    G4Hdf5AnalysisReader.hh
    G4Hdf5AnalysisReader.icc
    G4Hdf5FileManager.hh
    G4Hdf5HnFileManager.hh
    G4Hdf5HnFileManager.icc
    G4Hdf5HnRFileManager.hh
    G4Hdf5HnRFileManager.icc
    G4Hdf5NtupleManager.hh
    G4Hdf5NtupleFileManager.hh
    G4Hdf5RNtupleManager.hh
    G4Hdf5RFileManager.hh
    g4hdf5_defs.hh
  SOURCES
    G4Hdf5AnalysisManager.cc
    G4Hdf5AnalysisReader.cc
    G4Hdf5FileManager.cc
    G4Hdf5NtupleManager.cc
    G4Hdf5NtupleFileManager.cc
    G4Hdf5RFileManager.cc
    G4Hdf5RNtupleManager.cc)

geant4_module_link_libraries(G4hdf5 PUBLIC G4analysismng G4hntools G4globman G4tools ${HDF5_LIBRARIES})
