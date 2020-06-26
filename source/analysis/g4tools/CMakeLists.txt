#-----------------------------------------------------------------------
# Module : G4analysis
# Package: Geant4.src.G4analysis.g4tools
#
# CMakeLists.txt for g4tools, retained due to the nested headers. 
# The 'tools' subdirectory contains many headers that are needed and
# these are in a nested structure which needs to be replicated when we 
# install, i.e:
# <INCLUDE>
# +- geant4
#    +- tools
#       +- histo
#       |  +- h1d 
#       ...
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
  DESTINATION ${CMAKE_INSTALL_DATADIR}
  )
