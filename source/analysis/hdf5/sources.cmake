#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hdf5
# Package: Geant4.src.G4analysis.G4hdf5
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 27/07/2017
#
# $Id$
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${EXPAT_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/analysis/g4tools/include)
include_directories(${CMAKE_SOURCE_DIR}/source/analysis/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/analysis/hntools/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4xml
    HEADERS
        G4Hdf5AnalysisManager.hh
        G4Hdf5AnalysisManager.icc
        G4Hdf5AnalysisReader.hh
        G4Hdf5AnalysisReader.icc
        G4Hdf5FileManager.hh
        G4Hdf5NtupleManager.hh
        G4Hdf5RNtupleManager.hh
        G4Hdf5RFileManager.hh
        g4hdf5_defs.hh
        g4hdf5.hh
    SOURCES
        G4Hdf5AnalysisManager.cc
        G4Hdf5AnalysisReader.cc
        G4Hdf5FileManager.cc
        G4Hdf5NtupleManager.cc
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
