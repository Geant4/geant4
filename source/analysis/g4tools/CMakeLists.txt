#-----------------------------------------------------------------------
# CMakeLists.txt
# Module : G4analysis
# Package: Geant4.src.G4analysis.g4tools
#
# CMakeLists.txt for g4tools. We do not use the usual Geant4 sources.cmake
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
# $Id: CMakeLists.txt 72327 2013-07-16 06:52:43Z gcosmo $
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Ensure the header directory gets added to the list of ones to export
#
set_property(GLOBAL APPEND PROPERTY GEANT4_BUILDTREE_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  )

#-----------------------------------------------------------------------
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

