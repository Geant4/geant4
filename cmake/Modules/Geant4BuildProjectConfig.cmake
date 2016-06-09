# - Build Geant4Config.cmake file and support scripts for build and install.
#

#-----------------------------------------------------------------------
# Collect all global variables we need to export to the config files
# Do this here for now, later on we could collect them as we go.
#

# Compiler flags (because user apps are a bit dependent on them...)
set(GEANT4_COMPILER_FLAG_HINTS "#
set(Geant4_CXX_FLAGS \"${CMAKE_CXX_FLAGS}\")
set(Geant4_EXE_LINKER_FLAGS \"${CMAKE_EXE_LINKER_FLAGS}\")")

foreach(_mode DEBUG MINSIZEREL RELEASE RELWITHDEBINFO)
    set(GEANT4_COMPILER_FLAG_HINTS "${GEANT4_COMPILER_FLAG_HINTS}
set(Geant4_CXX_FLAGS_${_mode} \"${CMAKE_CXX_FLAGS_${_mode}}\")")
endforeach()


# Core compile definitions...
set(GEANT4_CORE_DEFINITIONS )

# Third party includes (libraries *should* be handled by the imports)
set(GEANT4_THIRD_PARTY_INCLUDES )

# Imports of third party packages used with imported targets
set(GEANT4_THIRD_PARTY_IMPORT_SETUP )

# Externals libraries that may be present
set(GEANT4_EXTERNALS_TARGETS )


# - Stuff from Geant4LibraryBuildOptions.cmake
if(GEANT4_BUILD_STORE_TRAJECTORY)
  list(APPEND GEANT4_CORE_DEFINITIONS -DG4_STORE_TRAJECTORY)
endif()

if(GEANT4_BUILD_VERBOSE_CODE)
  list(APPEND GEANT4_CORE_DEFINITIONS -DG4VERBOSE)
endif()

# - We do actually need G4LIB_BUILD_DLL on Windows, even for user 
# applications...
if(WIN32)
  list(APPEND GEANT4_CORE_DEFINITIONS -DG4LIB_BUILD_DLL)
endif()


# - Stuff from Geant4OptionalComponents.cmake
# - CLHEP
# If it's internal, add it to the externals list, if it's external, add the 
# include directories to the list of third party header paths
if(GEANT4_USE_SYSTEM_CLHEP)
  list(APPEND GEANT4_THIRD_PARTY_INCLUDES ${CLHEP_INCLUDE_DIRS})
else()
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

# - GDML
# Need to include Xerces-C headers becuase these do appear in the public
# interface of GDML. The library should then be in the LINK_INTERFACE of
# persistency...
if(GEANT4_USE_GDML)
  list(APPEND GEANT4_THIRD_PARTY_INCLUDES ${XERCESC_INCLUDE_DIRS})
endif()

# - Stuff from Geant4InterfaceOptions.cmake
if(GEANT4_USE_QT)
  list(APPEND GEANT4_THIRD_PARTY_INCLUDES 
    ${QT_INCLUDE_DIR}
    ${QT_QTCORE_INCLUDE_DIR}
    ${QT_QTGUI_INCLUDE_DIR}
    ${QT_QTOPENGL_INCLUDE_DIR}
    )

  # On WIN32, re-import the Qt targets.    
  if(WIN32)
    set(GEANT4_THIRD_PARTY_IMPORT_SETUP "${GEANT4_THIRD_PARTY_IMPORT_SETUP}
# Qt reimport on WIN32
set(QT_QMAKE_EXECUTABLE ${QT_QMAKE_EXECUTABLE})
set(QT_USE_IMPORTED_TARGETS ON)
find_package(Qt4 REQUIRED)"
        )
  endif()
endif()

#-----------------------------------------------------------------------
# - Generate Build Tree Configuration Files
#-----------------------------------------------------------------------
# Set needed variables for the build tree
set(GEANT4_CMAKE_DIR "${PROJECT_BINARY_DIR}")

# Set include path for build tree. This is always an absolute path, or 
# rather paths. We extract the paths from the global 
# GEANT4_BUILDTREE_INCLUDE_DIRS property and use this to create the 
# header setup
#
get_property(__geant4_buildtree_include_dirs GLOBAL PROPERTY
    GEANT4_BUILDTREE_INCLUDE_DIRS)

set(GEANT4_INCLUDE_DIR_SETUP "
# Geant4 configured for use from the build tree - absolute paths are used.
set(Geant4_INCLUDE_DIR ${__geant4_buildtree_include_dirs})
")

set(GEANT4_MODULE_PATH_SETUP "
# Geant4 configured for use CMake modules from source tree
set(CMAKE_MODULE_PATH \${CMAKE_MODULE_PATH} ${CMAKE_MODULE_PATH})
")

# Export targets from the build tree. We rely on the GEANT4_EXPORTED_TARGETS
# global property to list these for us.
#
get_property(__geant4_exported_targets GLOBAL PROPERTY GEANT4_EXPORTED_TARGETS)

export(TARGETS ${__geant4_exported_targets} FILE
    ${PROJECT_BINARY_DIR}/Geant4LibraryDepends.cmake)


# Configure the build tree config file...
configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/Geant4Config.cmake.in
    ${PROJECT_BINARY_DIR}/Geant4Config.cmake
    @ONLY)

# Configure the build tree versioning file
configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/Geant4ConfigVersion.cmake.in
    ${PROJECT_BINARY_DIR}/Geant4ConfigVersion.cmake
    @ONLY)

# Copy the custom modules into the build tree
configure_file(${PROJECT_SOURCE_DIR}/cmake/Modules/CMakeMacroParseArguments.cmake
  ${PROJECT_BINARY_DIR}/Modules/CMakeMacroParseArguments.cmake
  COPYONLY
  )

foreach(_mod AIDA HepMC Pythia6 ROOT StatTest)
  configure_file(${PROJECT_SOURCE_DIR}/cmake/Modules/Find${_mod}.cmake
    ${PROJECT_BINARY_DIR}/Modules/Find${_mod}.cmake
    COPYONLY
    )
endforeach()

# Copy the Main and Internal Use file into the build tree
configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/UseGeant4.cmake
    ${PROJECT_BINARY_DIR}/UseGeant4.cmake
    COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/UseGeant4_internal.cmake
    ${PROJECT_BINARY_DIR}/UseGeant4_internal.cmake
    COPYONLY)


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
get_filename_component(Geant4_INCLUDE_DIR \"\${_thisdir}/${GEANT4_RELATIVE_HEADER_PATH}\" ABSOLUTE)
")
endif()

set(GEANT4_MODULE_PATH_SETUP)

# Install exported targets file for the install tree - we just install
# the named export
install(EXPORT Geant4LibraryDepends 
    DESTINATION ${GEANT4_CMAKE_DIR}
    COMPONENT Development
)


# Configure the install tree config file...
configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/Geant4Config.cmake.in
    ${PROJECT_BINARY_DIR}/InstallTreeFiles/Geant4Config.cmake
    @ONLY
)

# Configure the install tree config versioning file...
configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/Geant4ConfigVersion.cmake.in
    ${PROJECT_BINARY_DIR}/InstallTreeFiles/Geant4ConfigVersion.cmake
    @ONLY
)

# Install the config, config versioning and use files
install(FILES
  ${PROJECT_BINARY_DIR}/InstallTreeFiles/Geant4Config.cmake
  ${PROJECT_BINARY_DIR}/InstallTreeFiles/Geant4ConfigVersion.cmake
  ${PROJECT_SOURCE_DIR}/cmake/Templates/UseGeant4.cmake
  DESTINATION ${GEANT4_CMAKE_DIR}
  COMPONENT Development
  )

# Install the custom modules for the examples
install(DIRECTORY
  ${PROJECT_BINARY_DIR}/Modules
  DESTINATION ${GEANT4_CMAKE_DIR}
  COMPONENT Development
  )

