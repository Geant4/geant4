#------------------------------------------------------------------------------
# sources.cmake
# Module : G4csv
# Package: Geant4.src.G4analysis.G4xml
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 15/07/2013
#
# $Id: sources.cmake 105338 2017-07-21 09:14:27Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${EXPAT_INCLUDE_DIRS})

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
        G4XmlAnalysisManager.hh
        G4XmlAnalysisManager.icc
        G4XmlAnalysisReader.hh
        G4XmlAnalysisReader.icc
        G4XmlFileManager.hh
        G4XmlNtupleManager.hh
        G4XmlRFileManager.hh
        G4XmlRNtupleManager.hh
        g4xml_defs.hh
        g4xml.hh
    SOURCES
        G4XmlAnalysisManager.cc
        G4XmlAnalysisReader.cc
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
