#------------------------------------------------------------------------------
# sources.cmake
# Module : G4csv
# Package: Geant4.src.G4analysis.G4analysismng
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 15/07/2013
#
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${ZLIB_INCLUDE_DIRS})
if (GEANT4_USE_HDF5)
  include_directories(${HDF5_INCLUDE_DIRS})
  set(G4analysisfac_LINK_LIBRARIES Geant4::HDF5)
endif()
if(GEANT4_USE_FREETYPE)
  set(G4analysisfac_LINK_LIBRARIES Freetype::Freetype ${G4analysisfac_LINK_LIBRARIES})
endif()

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/analysis/g4tools/include)
include_directories(${CMAKE_SOURCE_DIR}/source/analysis/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/analysis/hntools/include)
include_directories(${CMAKE_SOURCE_DIR}/source/analysis/csv/include)
include_directories(${CMAKE_SOURCE_DIR}/source/analysis/root/include)
include_directories(${CMAKE_SOURCE_DIR}/source/analysis/xml/include)
if (GEANT4_USE_HDF5)
include_directories(${CMAKE_SOURCE_DIR}/source/analysis/hdf5/include)
endif()

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4analysisfac
    HEADERS
        g4analysis.hh
        g4analysis_defs.hh
    SOURCES
        g4analysis.cc
    GRANULAR_DEPENDENCIES
        G4globman
        G4intercoms
    GLOBAL_DEPENDENCIES
        G4global
        G4intercoms
    LINK_LIBRARIES
        ${G4analysisfac_LINK_LIBRARIES}
)

# List any source specific properties here
