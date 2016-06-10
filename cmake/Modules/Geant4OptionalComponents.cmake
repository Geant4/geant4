# - Setup core required and optional components of Geant4
#
# Here we provide options to enable and configure required and optional
# components, which may require third party libraries.
#
# We don't configure User Interface options here because these require
# a higher degree of configuration so to keep things neat these have their
# own Module.
#
# Options configured here:
#  CLHEP   - Control use of internal G4clhep, or locate external CLHEP
#            Also control selection of singular or modular CLHEP libs
#  EXPAT   - Control use of internal G4expat, or locate external EXPAT.
#  ZLIB    - Control use of internal G4zlib, or locate external ZLIB
#  GDML    - Requires external XercesC
#  G3TOG4  - UNIX only
#  USOLIDS - Allow use of USolids classes in geometry, using
#            internal Usolids by default, plus option to use system
#            version

#-----------------------------------------------------------------------
# Find required CLHEP package
# We prefer to use our internal CLHEP library, but provide an option
# to allow an external CLHEP to be used should the user require this.
# We also allow that it can be automatically enabled by providing
# the CLHEP_ROOT_DIR option (which FindCLHEP will recognize)
#
# As requested by ATLAS, an additional option for preferring use of
# CLHEP's granular libs is provided when using a system CLHEP.
#
# KNOWNISSUE : For internal CLHEP, how to deal with static and shared?
if(CLHEP_ROOT_DIR)
  set(_default_use_system_clhep ON)
else()
  set(_default_use_system_clhep OFF)
endif()

option(GEANT4_USE_SYSTEM_CLHEP "Use system CLHEP library" ${_default_use_system_clhep})

if(GEANT4_USE_SYSTEM_CLHEP)
  set(__system_clhep_mode " (singular)")
  # Further advanced option to select granular CLHEP libs
  option(GEANT4_USE_SYSTEM_CLHEP_GRANULAR "Use system CLHEP granular libraries" OFF)
  mark_as_advanced(GEANT4_USE_SYSTEM_CLHEP_GRANULAR)

  if(GEANT4_USE_SYSTEM_CLHEP_GRANULAR)
    set(__g4_clhep_components
      Evaluator
      Geometry
      Random
      Vector
      )
    set(__system_clhep_mode " (granular)")
  endif()

  # We keep this as required, because if the user chooses to use a
  # system option we assume that we absolutely, positively require this.
  find_package(CLHEP 2.1.2.3 REQUIRED ${__g4_clhep_components})
  set(GEANT4_USE_SYSTEM_CLHEP TRUE)
else()
  set(CLHEP_FOUND TRUE)
  set(GEANT4_USE_BUILTIN_CLHEP TRUE)
  set(CLHEP_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/source/externals/clhep/include")
  if(BUILD_SHARED_LIBS)
    set(CLHEP_LIBRARIES G4clhep)
  else()
    set(CLHEP_LIBRARIES G4clhep-static)
  endif()
endif()

GEANT4_ADD_FEATURE(GEANT4_USE_SYSTEM_CLHEP "Using system CLHEP library${__system_clhep_mode}")

#-----------------------------------------------------------------------
# Find required EXPAT package
# We always use the internal G4expat on WIN32.
# On other platforms, we default to use the system library.
# If we use the internal G4expat, fix the EXPAT_XXX variables to point
# to the internal location of expat headers and library target so Geant4
# clients of expat have a consistent interface.
#
# KNOWNISSUE : For internal expat, how to deal with static and shared?
if(WIN32)
  set(EXPAT_FOUND TRUE)
  set(GEANT4_USE_BUILTIN_EXPAT TRUE)
  set(EXPAT_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/source/externals/expat/include)
  set(EXPAT_LIBRARIES G4expat)
else()
  option(GEANT4_USE_SYSTEM_EXPAT "Use system Expat library" ON)

  if(GEANT4_USE_SYSTEM_EXPAT)
    # If system package requested, make damn sure we find it
    find_package(EXPAT REQUIRED)
  else()
    set(EXPAT_FOUND TRUE)
    set(GEANT4_USE_BUILTIN_EXPAT TRUE)
    set(EXPAT_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/source/externals/expat/include)
    if(BUILD_SHARED_LIBS)
      set(EXPAT_LIBRARIES G4expat)
    else()
      set(EXPAT_LIBRARIES G4expat-static)
    endif()
  endif()
endif()

GEANT4_ADD_FEATURE(GEANT4_USE_SYSTEM_EXPAT "Using system EXPAT library")

#-----------------------------------------------------------------------
# Find required ZLIB package
# Default to use internal zlib, otherwise point interface variables to
# internal zlib
option(GEANT4_USE_SYSTEM_ZLIB "Use system zlib library" OFF)
if(GEANT4_USE_SYSTEM_ZLIB)
  find_package(ZLIB REQUIRED)

  # NB : FindZLIB on cmake < 2.8 does not set the ZLIB_INCLUDE_DIRS
  # variable, only the ZLIB_INCLUDE_DIR variable. Set the DIRS variable
  # here for backward compatibility.
  if(${CMAKE_VERSION} VERSION_LESS "2.8.0")
    set(ZLIB_INCLUDE_DIRS "${ZLIB_INCLUDE_DIR}")
  endif()
else()
  set(ZLIB_FOUND TRUE)
  set(GEANT4_USE_BUILTIN_ZLIB TRUE)
  set(ZLIB_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/source/externals/zlib/include
                        ${PROJECT_BINARY_DIR}/source/externals/zlib)
  if(BUILD_SHARED_LIBS)
    set(ZLIB_LIBRARIES G4zlib)
  else()
    set(ZLIB_LIBRARIES G4zlib-static)
  endif()
endif()

GEANT4_ADD_FEATURE(GEANT4_USE_SYSTEM_ZLIB "Using system zlib library")

#-----------------------------------------------------------------------
# Optional Support for GDML - requires Xerces-C package
#
if(XERCESC_ROOT_DIR)
  set(_default_use_gdml ON)
else()
  set(_default_use_gdml OFF)
endif()

option(GEANT4_USE_GDML "Build Geant4 with GDML support" ${_default_use_gdml}
)

if(GEANT4_USE_GDML)
  find_package(XercesC REQUIRED)
endif()

GEANT4_ADD_FEATURE(GEANT4_USE_GDML "Building Geant4 with GDML support")

#-----------------------------------------------------------------------
# Optional support for G3TOG4 convertion interface.
# We do not build the rztog4 application.
# -- OLDER NOTES --
# The G3toG4 *library* should always be built, but the rztog4 application
# requires a Fortran compiler AND CERNLIB, so is optional.
# Only on *NIX because converter requires CERNLIB, and Windows support for
# this is an unknown quantity at present (can always change later).
# -- OLDER NOTES --
#
if(UNIX)
  option(GEANT4_USE_G3TOG4 "Build Geant3 ASCII call list reader library" OFF)
  GEANT4_ADD_FEATURE(GEANT4_USE_G3TOG4 "Building Geant3 ASCII call list reader library")
endif()

#-----------------------------------------------------------------------
# Optional replacement
# - Advanced option only
# - If enabled, require
#   1) Global Compile definition G4GEOM_USE_USOLIDS (also exported)
#   2) Use internal USolids by default, otherwise searching for system
#      install
option(GEANT4_USE_USOLIDS "EXPERIMENTAL: Replace Geant4 solids with USolids equivalents" OFF)
option(GEANT4_USE_SYSTEM_USOLIDS "Use system USolids library" OFF)
mark_as_advanced(GEANT4_USE_USOLIDS GEANT4_USE_SYSTEM_USOLIDS)

# - G4USolids setup
if(GEANT4_USE_USOLIDS)
  add_definitions(-DG4GEOM_USE_USOLIDS)
endif()

# - Internal or external USolids
if(GEANT4_USE_SYSTEM_USOLIDS)
  find_package(USolids REQUIRED)
else()
  set(USolids_FOUND TRUE)
  set(GEANT4_USE_BUILTIN_USOLIDS TRUE)
  set(USOLIDS_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/source/externals/usolids/include)

  if(BUILD_SHARED_LIBS)
    set(USOLIDS_LIBRARIES G4geomUSolids)
  else()
    set(USOLIDS_LIBRARIES G4geomUSolids-static)
  endif()
  # Include dirs here because of the large number of G4 users
  # of solids
  # Resolve once we use Modern INTERFACE_INCLUDES
  include_directories(${USOLIDS_INCLUDE_DIRS})
endif()

GEANT4_ADD_FEATURE(GEANT4_USE_USOLIDS "Replacing Geant4 solids with USolids equivalents (EXPERIMENTAL)")
GEANT4_ADD_FEATURE(GEANT4_USE_SYSTEM_USOLIDS "Using system USolids library")

