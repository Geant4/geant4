#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4analysis
# Package: Geant4.src.G4analysis
#
# CMakeLists.txt for G4analysis. We do not use the usual Geant4 sources.cmake
# approach because G4analysis requires a custom install solution for its
# headers. The 'tools' subdirectory contains many headers that are needed and
# these are in a nested structure which needs to be replicated when we 
# install, i.e:
# <INCLUDE>
# +- geant4
#    +- tools
#       +- histo
#       |  +- h1d 
#       ...
#
# Nevertheless, we're able to build the library via the standard
# GEANT4_LIBRARY_TARGET macro.
#
# Created on: 28/07/2011
#
# $Date$
# $Revision$
# $Author$
#
#------------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Includes for this modules
#
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)

# List the headers and sources
#
# - Headers
# Need to add tools, but since these don't use a standard naming scheme,
# will need to work out how to either locate them, or generate a list.
# This isn't critical for building the analysis library.
set(G4analysis_HEADERS
    include/G4AnalysisMessenger.hh
    include/G4AnalysisVerbose.hh
    include/G4CsvAnalysisManager.hh
    include/g4csv_defs.hh
    include/g4csv.hh
    include/G4Fcn.hh
    include/G4HnInformation.hh
    include/G4RootAnalysisManager.hh
    include/g4root_defs.hh
    include/g4root.hh
    include/G4VAnalysisManager.hh
    include/G4XmlAnalysisManager.hh
    include/g4xml_defs.hh
    include/g4xml.hh
)

# - Sources
set(G4analysis_SOURCES
    src/G4AnalysisMessenger.cc
    src/G4AnalysisVerbose.cc
    src/G4CsvAnalysisManager.cc
    src/G4RootAnalysisManager.cc
    src/G4VAnalysisManager.cc
    src/G4XmlAnalysisManager.cc
)


#----------------------------------------------------------------------------
# Now add the library targets
#
include(Geant4MacroLibraryTargets)
if(GEANT4_BUILD_GRANULAR_LIBS)
  GEANT4_LIBRARY_TARGET(NAME G4analysis 
    SOURCES ${G4analysis_SOURCES} ${G4analysis_HEADERS}
    GEANT4_LINK_LIBRARIES G4globman G4intercoms
    )
else()
  GEANT4_LIBRARY_TARGET(NAME G4analysis 
    SOURCES ${G4analysis_SOURCES} ${G4analysis_HEADERS}
    GEANT4_LINK_LIBRARIES G4global G4intercoms
    )
endif()

#----------------------------------------------------------------------------
# Ensure the header directory gets added to the list of ones to export
#
set_property(GLOBAL APPEND PROPERTY GEANT4_BUILDTREE_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)


#----------------------------------------------------------------------------
# Install the headers, retaining their nested structure
#
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    COMPONENT Development
    PATTERN ".svn" EXCLUDE
)

#----------------------------------------------------------------------------
# We also need to install the tools license
#
install(FILES tools.license
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Geant4-${Geant4_VERSION}
)

