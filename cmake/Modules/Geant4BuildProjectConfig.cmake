# - Build Geant4Config.cmake file and support scripts for build and install.
#

#----------------------------------------------------------------------------
# - Generate Build Tree Configuration Files
#
# Export targets from the build tree. We rely on the GEANT4_EXPORTED_TARGETS
# global property to list these for us.
#
get_property(__geant4_exported_targets GLOBAL PROPERTY GEANT4_EXPORTED_TARGETS)

export(TARGETS ${__geant4_exported_targets} FILE
    ${PROJECT_BINARY_DIR}/Geant4LibraryDepends.cmake)

# Set include path for build tree. This is always an absolute path, or rather
# paths. We extract the paths from the global GEANT4_BUILDTREE_INCLUDE_DIRS
# property and use this to create the header setup
#
get_property(__geant4_buildtree_include_dirs GLOBAL PROPERTY
    GEANT4_BUILDTREE_INCLUDE_DIRS)

set(GEANT4_INCLUDE_DIR_SETUP "
# Geant4 configured for use from the build tree - absolute paths are used.
set(Geant4_INCLUDE_DIR ${__geant4_buildtree_include_dirs})
")

# Configure the build tree config file...
configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/Geant4Config.cmake.in
    ${PROJECT_BINARY_DIR}/Geant4Config.cmake
    @ONLY)

# Configure the build tree versioning file
configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/Geant4ConfigVersion.cmake.in
    ${PROJECT_BINARY_DIR}/Geant4ConfigVersion.cmake
    @ONLY)




#----------------------------------------------------------------------------
# - Generate Install Tree Configuration Files
#

#----------------------------------------------------------------------------
# Install exported targets file for the install tree - we just install
# the named export
#
install(EXPORT Geant4LibraryDepends
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}-${${PROJECT_NAME}_VERSION}
    COMPONENT Development)





