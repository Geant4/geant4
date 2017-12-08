#------------------------------------------------------------------------------
# sources.cmake
# Module : G4csv
# Package: Geant4.src.G4analysis.G4csv
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
GEANT4_DEFINE_MODULE(NAME G4csv
    HEADERS
        G4CsvAnalysisManager.hh
        G4CsvAnalysisManager.icc
        G4CsvAnalysisReader.hh
        G4CsvAnalysisReader.icc
        G4CsvFileManager.hh
        G4CsvNtupleManager.hh
        G4CsvRFileManager.hh
        G4CsvRNtupleManager.hh
        g4csv_defs.hh
        g4csv.hh
    SOURCES
        G4CsvAnalysisManager.cc
        G4CsvAnalysisReader.cc
        G4CsvFileManager.cc
        G4CsvNtupleManager.cc
        G4CsvRFileManager.cc
        G4CsvRNtupleManager.cc
    GRANULAR_DEPENDENCIES
        G4globman
        G4intercoms
        G4analysismng
        G4hntools
    GLOBAL_DEPENDENCIES
        G4global
        G4intercoms
    LINK_LIBRARIES
)

# List any source specific properties here
