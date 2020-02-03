#.rst:
# G4ConfigureCMakeHelpers
# -----------------------
#
# This module configures and installs CMake modules allowing clients
# to find and use Geant4 libraries using CMake's find_package command.
#

#-----------------------------------------------------------------
# License and Disclaimer
#
# The  Geant4 software  is  copyright of the Copyright Holders  of
# the Geant4 Collaboration.  It is provided  under  the terms  and
# conditions of the Geant4 Software License,  included in the file
# LICENSE and available at  http://cern.ch/geant4/license .  These
# include a list of copyright holders.
#
# Neither the authors of this software system, nor their employing
# institutes,nor the agencies providing financial support for this
# work  make  any representation or  warranty, express or implied,
# regarding  this  software system or assume any liability for its
# use.  Please see the license in the file  LICENSE  and URL above
# for the full disclaimer and the limitation of liability.
#
# This  code  implementation is the result of  the  scientific and
# technical work of the GEANT4 collaboration.
# By using,  copying,  modifying or  distributing the software (or
# any work based  on the software)  you  agree  to acknowledge its
# use  in  resulting  scientific  publications,  and indicate your
# acceptance of all terms of the Geant4 Software license.
#
#-----------------------------------------------------------------

#-----------------------------------------------------------------------
# Collect all global variables we need to export to the config files
# Do this here for now, later on we could collect them as we go.
#

# Compiler flags (because user apps are a bit dependent on them...)
set(GEANT4_COMPILER_FLAG_HINTS "#
set(Geant4_CXX_FLAGS \"${CMAKE_CXX_FLAGS} ${GEANT4_CXXSTD_FLAGS}\")
set(Geant4_EXE_LINKER_FLAGS \"${CMAKE_EXE_LINKER_FLAGS}\")")

foreach(_mode DEBUG MINSIZEREL RELEASE RELWITHDEBINFO)
  set(GEANT4_COMPILER_FLAG_HINTS "${GEANT4_COMPILER_FLAG_HINTS}
set(Geant4_CXX_FLAGS_${_mode} \"${CMAKE_CXX_FLAGS_${_mode}}\")")
endforeach()

if(NOT CMAKE_CONFIGURATION_TYPES)
  set(GEANT4_COMPILER_FLAG_HINTS "${GEANT4_COMPILER_FLAG_HINTS}
set(Geant4_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")")
endif()

# Core compile definitions...
set(GEANT4_CORE_DEFINITIONS )

# Third party includes (libraries *should* be handled by the imports)
set(GEANT4_THIRD_PARTY_INCLUDES )

# Imports of third party packages used with imported targets
set(GEANT4_THIRD_PARTY_IMPORT_SETUP )

# Externals libraries that may be present
set(GEANT4_EXTERNALS_TARGETS )

# - Stuff from Geant4OptionalComponents.cmake
# - CLHEP
# If it's internal, add it to the externals list
if(NOT GEANT4_USE_SYSTEM_CLHEP)
  list(APPEND GEANT4_EXTERNALS_TARGETS G4clhep)
endif()

# - Expat
# If it's internal, add it to the externals list
if(NOT GEANT4_USE_SYSTEM_EXPAT)
  list(APPEND GEANT4_EXTERNALS_TARGETS G4expat)
endif()

# - ZLIB
# If it's internal, add it to the externals list
if(NOT GEANT4_USE_SYSTEM_ZLIB)
  list(APPEND GEANT4_EXTERNALS_TARGETS G4zlib)
endif()

# - USolids
# Compile definitions
if(GEANT4_USE_USOLIDS OR GEANT4_USE_PARTIAL_USOLIDS)
  set(GEANT4_USE_USOLIDS_EITHER ON)
endif()

# - Stuff from Geant4InterfaceOptions.cmake

#-----------------------------------------------------------------------
# - Common Build/Install Tree Configuration files
#-----------------------------------------------------------------------
# External package variables
geant4_export_package_variables(${PROJECT_BINARY_DIR}/Geant4PackageCache.cmake)

# Versioning file
configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/Templates/Geant4ConfigVersion.cmake.in
  ${PROJECT_BINARY_DIR}/Geant4ConfigVersion.cmake
  @ONLY
  )

#-----------------------------------------------------------------------
# - Generate Build Tree Configuration Files
#-----------------------------------------------------------------------
# Set needed variables for the build tree
set(GEANT4_CMAKE_DIR "${PROJECT_BINARY_DIR}")

# Set include path for build tree. This is always an absolute path, or
# rather paths. We extract the paths from the global
# GEANT4_BUILDTREE_INCLUDE_DIRS property and use this to create the
# header setup
# This is *ONLY* required in case clients use ROOT dictionaries, whose
# generation mechanism is blind to usage requirements
#
get_property(__geant4_buildtree_include_dirs GLOBAL PROPERTY
  GEANT4_BUILDTREE_INCLUDE_DIRS
  )

set(GEANT4_INCLUDE_DIR_SETUP "
# Geant4 configured for use from the build tree - absolute paths are used.
set(Geant4_INCLUDE_DIR \"${__geant4_buildtree_include_dirs}\")
")

# Geant4 data used in build tree
geant4_export_datasets(BUILD GEANT4_DATASET_DESCRIPTIONS)

# Export targets in the Geant4LibraryDepends export set from the build tree
export(EXPORT Geant4LibraryDepends
  NAMESPACE Geant4::
  FILE ${PROJECT_BINARY_DIR}/Geant4LibraryDepends.cmake
  )

# Configure the build tree config file...
configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/Templates/Geant4Config.cmake.in
  ${PROJECT_BINARY_DIR}/Geant4Config.cmake
  @ONLY
  )

# Copy the custom modules into the build tree
configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/Modules/IntelCompileFeatures.cmake
  ${PROJECT_BINARY_DIR}/Modules/IntelCompileFeatures.cmake
  COPYONLY
  )

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/Modules/MSVCCompileFeatures.cmake
  ${PROJECT_BINARY_DIR}/Modules/MSVCCompileFeatures.cmake
  COPYONLY
  )

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/Modules/G4EXPATShim.cmake
  ${PROJECT_BINARY_DIR}/G4EXPATShim.cmake
  COPYONLY
  )

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/Modules/G4FreetypeShim.cmake
  ${PROJECT_BINARY_DIR}/G4FreetypeShim.cmake
  COPYONLY
  )

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/Modules/G4HDF5Shim.cmake
  ${PROJECT_BINARY_DIR}/G4HDF5Shim.cmake
  COPYONLY
  )

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/Modules/G4MotifShim.cmake
  ${PROJECT_BINARY_DIR}/G4MotifShim.cmake
  COPYONLY
)

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/Modules/G4VecGeomShim.cmake
  ${PROJECT_BINARY_DIR}/G4VecGeomShim.cmake
  COPYONLY
  )

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/Modules/G4X11Shim.cmake
  ${PROJECT_BINARY_DIR}/G4X11Shim.cmake
  COPYONLY
)


foreach(_mod AIDA HepMC Pythia6 StatTest TBB XQuartzGL)
  configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/Modules/Find${_mod}.cmake
    ${PROJECT_BINARY_DIR}/Modules/Find${_mod}.cmake
    COPYONLY
    )
endforeach()

# Copy the Main and Internal Use file into the build tree
configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/Templates/UseGeant4.cmake
  ${PROJECT_BINARY_DIR}/UseGeant4.cmake
  COPYONLY
  )

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/Templates/UseGeant4_internal.cmake
  ${PROJECT_BINARY_DIR}/UseGeant4_internal.cmake
  COPYONLY
  )

#-----------------------------------------------------------------------
# - Generate Install Tree Configuration Files
#-----------------------------------------------------------------------
# Set needed variables for the install tree
set(GEANT4_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}-${${PROJECT_NAME}_VERSION})

# Header path for install tree is dependent on whether we have a relocatable
# install.
if(CMAKE_INSTALL_IS_NONRELOCATABLE)
  # Use ABSOLUTE paths...
  set(GEANT4_INCLUDE_DIR_SETUP "
# Geant4 configured for the install tree with absolute paths, so use these
set(Geant4_INCLUDE_DIR \"${CMAKE_INSTALL_FULL_INCLUDEDIR}/Geant4\")
  ")
else()
  # Use RELATIVE paths... Where we measure relative to GEANT4_CMAKE_DIR
  file(RELATIVE_PATH GEANT4_RELATIVE_HEADER_PATH
    ${CMAKE_INSTALL_PREFIX}/${GEANT4_CMAKE_DIR}
    ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    )

  set(GEANT4_INCLUDE_DIR_SETUP "
# Geant4 configured for the install with relative paths, so use these
get_filename_component(Geant4_INCLUDE_DIR \"\${_geant4_thisdir}/${GEANT4_RELATIVE_HEADER_PATH}\" ABSOLUTE)
  ")
endif()

# Geant4 data used in install tree
geant4_export_datasets(INSTALL GEANT4_DATASET_DESCRIPTIONS)

# Install exported targets file for the install tree - we just install
# the named export
install(EXPORT Geant4LibraryDepends
  NAMESPACE Geant4::
  DESTINATION ${GEANT4_CMAKE_DIR}
  COMPONENT Development
  )

# Configure the install tree config file...
configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/Templates/Geant4Config.cmake.in
  ${PROJECT_BINARY_DIR}/InstallTreeFiles/Geant4Config.cmake
  @ONLY
  )

# Install the config, config versioning and use files
install(FILES
  ${PROJECT_BINARY_DIR}/InstallTreeFiles/Geant4Config.cmake
  ${PROJECT_BINARY_DIR}/Geant4ConfigVersion.cmake
  ${PROJECT_BINARY_DIR}/G4EXPATShim.cmake
  ${PROJECT_BINARY_DIR}/G4FreetypeShim.cmake
  ${PROJECT_BINARY_DIR}/G4HDF5Shim.cmake
  ${PROJECT_BINARY_DIR}/G4VecGeomShim.cmake
  ${PROJECT_BINARY_DIR}/G4MotifShim.cmake
  ${PROJECT_BINARY_DIR}/G4X11Shim.cmake
  ${PROJECT_SOURCE_DIR}/cmake/Templates/UseGeant4.cmake
  DESTINATION ${GEANT4_CMAKE_DIR}
  COMPONENT Development
  )

# Install the package settings file if required (always for now)
option(GEANT4_INSTALL_PACKAGE_CACHE "Install file recording build-time locations of required packages" ON)
mark_as_advanced(GEANT4_INSTALL_PACKAGE_CACHE)
if(GEANT4_INSTALL_PACKAGE_CACHE)
  install(FILES ${PROJECT_BINARY_DIR}/Geant4PackageCache.cmake
    DESTINATION ${GEANT4_CMAKE_DIR}
    COMPONENT Development
    )
endif()

# Install the custom modules for dependencies/examples
install(DIRECTORY
  ${PROJECT_BINARY_DIR}/Modules
  DESTINATION ${GEANT4_CMAKE_DIR}
  COMPONENT Development
  )

