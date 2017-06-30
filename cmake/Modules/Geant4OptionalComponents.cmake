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
#  CLHEP    - Control use of internal G4clhep, or locate external CLHEP
#             Also control selection of singular or modular CLHEP libs
#  EXPAT    - Control use of internal G4expat, or locate external EXPAT.
#  ZLIB     - Control use of internal G4zlib, or locate external ZLIB
#  GDML     - Requires external XercesC
#  G3TOG4   - UNIX only
#  USOLIDS  - Allow use of USolids classes in geometry, and must find
#             external USolids install
#  FREETYPE - For analysis module, find external FreeType install

#-----------------------------------------------------------------------
# CLHEP
# ^^^^^
#
# By default, Geant4 is built with its own internal copy of the CLHEP
# libraries and headers. An external install of CLHEP may be used
# instead by setting the `GEANT4_USE_SYSTEM_CLHEP` to `ON`.
#
# CMake can be pointed to the required external install of CLHEP
# by providing the `CLHEP_ROOT_DIR` option. Setting this variable
# will automatically enable use of an external CLHEP.
#
# When using an external install of CLHEP, the default is to link
# to the full CLHEP library. An additional advanced option
# `GEANT4_USE_SYSTEM_CLHEP_GRANULAR` is available to configure
# use of the modular CLHEP libraries.
#
if(CLHEP_ROOT_DIR)
  set(_default_use_system_clhep ON)
  # CLHEP_ROOT_DIR is non-standard for Config mode, so attempt temporary translation
  # into CMAKE_PREFIX_PATH. Could also add as PATHS to find_package,
  # but CMAKE_PREFIX_PATH has higher search priority.
  list(INSERT CMAKE_PREFIX_PATH 0 "${CLHEP_ROOT_DIR}")
else()
  set(_default_use_system_clhep OFF)
endif()

option(GEANT4_USE_SYSTEM_CLHEP "Use system CLHEP library" ${_default_use_system_clhep})
cmake_dependent_option(GEANT4_USE_SYSTEM_CLHEP_GRANULAR
  "Use system CLHEP granular libraries" OFF
  "GEANT4_USE_SYSTEM_CLHEP" OFF
  )
mark_as_advanced(GEANT4_USE_SYSTEM_CLHEP_GRANULAR)

if(GEANT4_USE_SYSTEM_CLHEP)
  set(__system_clhep_mode " (singular)")

  if(GEANT4_USE_SYSTEM_CLHEP_GRANULAR)
    set(__g4_clhep_components
      Evaluator
      Geometry
      Random
      Vector
      )
    set(__system_clhep_mode " (granular)")
  endif()

  # Find CLHEP using config-mode (i.e. CLHEPConfig)
  # Note that the imported targets have different properties
  # depending on version:
  # >= 2.3.3.0 Imported targets have INTERFACE_INCLUDE_DIRECTORIES
  # >= 2.3.1.0 Imported targets are namespaced, all components available
  # < 2.3.1.0 Only CLHEP/CLHEPS available as imported targets
  #
  # By supplying components, CLHEP_LIBRARIES is populated with the
  # appropriate set of imported targets.
  find_package(CLHEP 2.3.3.0 REQUIRED ${__g4_clhep_components} CONFIG)
else()
  set(CLHEP_FOUND TRUE)
  # TODO: CLHEP_INCLUDE_DIRS still required for windows support
  #       Current way of building DLLs requires build of a temporary
  #       archive lib, and this does not link to dependencies, hence
  #       it does not get their usage requirements.
  #       May be fixable by linking or use of better DLL build mechanism
  #       available in CMake >= 3.4
  #       As other externals will need the same treatment, return
  #       to old CLHEP_INCLUDE_DIRS setting
  #
  set(CLHEP_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/source/externals/clhep/include")

  # All G4*** targets are handled internally to use static or shared
  # as appropriate, so only need to declare core target name
  set(CLHEP_LIBRARIES G4clhep)
endif()

geant4_add_feature(GEANT4_USE_SYSTEM_CLHEP "Using system CLHEP library${__system_clhep_mode}")

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
# Optional support for USolids - Requires USolids package
# GEANT4_USE_USOLIDS is a list argument which must take values:
# - Empty/CMake boolean false : DEFAULT, do not replace G4 solids with USolids
# - ALL/CMake boolean true : Replace all G4 solids with USolids
# - List : Replace listed G4 solids with Usolids equivalents. Valid elements are
set(GEANT4_USOLIDS_SHAPES
  BOX
  CONS
  EXTRUDEDSOLID
  GENERICPOLYCONE
  GENERICTRAP
  ORB
  PARABOLOID
  POLYCONE
  POLYHEDRA
  SPHERE
  TET
  TRAP
  TRD
  TORUS
  TUBS
  )

set(GEANT4_USE_USOLIDS OFF CACHE STRING "EXPERIMENTAL: List Geant4 solids to replace with USolids equivalents (ALL;${GEANT4_USOLIDS_SHAPES})")
mark_as_advanced(GEANT4_USE_USOLIDS)

# - Must be a CMake list - no real way to enforce this
# Preprocess by stripping any spaces and uppercasing
string(REPLACE " " "" __g4_usolids_shape_list "${GEANT4_USE_USOLIDS}")
string(TOUPPER "${__g4_usolids_shape_list}" __g4_usolids_shape_list)

# - Check for use of all or usolids
set(GEANT4_USE_ALL_USOLIDS OFF)
string(REGEX MATCH "ALL|^(1|ON|YES|TRUE)" GEANT4_USE_ALL_USOLIDS "${__g4_usolids_shape_list}")

# - Check and validate partial non-empty list
if(NOT GEANT4_USE_ALL_USOLIDS AND GEANT4_USE_USOLIDS)
  set(GEANT4_USE_PARTIAL_USOLIDS ON)
  set(GEANT4_USE_PARTIAL_USOLIDS_SHAPE_LIST)

  foreach(__g4_usolids_requested_shape ${__g4_usolids_shape_list})
    list(FIND GEANT4_USOLIDS_SHAPES "${__g4_usolids_requested_shape}" __g4solids_shape_index)
    if(__g4solids_shape_index GREATER -1)
      list(APPEND GEANT4_USE_PARTIAL_USOLIDS_SHAPE_LIST ${__g4_usolids_requested_shape})
    else()
      message(FATAL_ERROR "GEANT4_USE_USOLIDS: unknown shape '${__g4_usolids_requested_shape}' in input, must be one of\n${GEANT4_USOLIDS_SHAPES}\n")
    endif()
  endforeach()
endif()


# - G4USolids setup
if(GEANT4_USE_ALL_USOLIDS OR GEANT4_USE_PARTIAL_USOLIDS)
  find_package(USolids REQUIRED)

  if(GEANT4_USE_ALL_USOLIDS)
    set(GEANT4_USOLIDS_COMPILE_DEFINITIONS "-DG4GEOM_USE_USOLIDS")
    GEANT4_ADD_FEATURE(GEANT4_USE_USOLIDS "Replacing all Geant4 solids with USolids equivalents (EXPERIMENTAL)")
  else()
    set(GEANT4_USOLIDS_COMPILE_DEFINITIONS "-DG4GEOM_USE_PARTIAL_USOLIDS")
    foreach(__g4_usolid_shape ${GEANT4_USE_PARTIAL_USOLIDS_SHAPE_LIST})
      list(APPEND GEANT4_USOLIDS_COMPILE_DEFINITIONS "-DG4GEOM_USE_U${__g4_usolid_shape}")
    endforeach()
    GEANT4_ADD_FEATURE(GEANT4_USE_USOLIDS "Replacing Geant4 solids with USolids equivalents for ${GEANT4_USE_PARTIAL_USOLIDS_SHAPE_LIST} (EXPERIMENTAL)")
  endif()
  list (APPEND GEANT4_USOLIDS_COMPILE_DEFINITIONS ${VECGEOM_DEFINITIONS})

  # Combined definitions
  add_definitions(${GEANT4_USOLIDS_COMPILE_DEFINITIONS})

  # Add USolids inc dirs here - can be removed once USolids supports
  # INTERFACE_INCLUDE_DIRECTORIES
  include_directories(${USOLIDS_INCLUDE_DIRS} ${VECGEOM_EXTERNAL_INCLUDES})
endif()


#-----------------------------------------------------------------------
# Optional support for Freetype - Requires external Freetype install
#
option(GEANT4_USE_FREETYPE "Build Geant4 analysis library with Freetype support" OFF)
mark_as_advanced(GEANT4_USE_FREETYPE)

if(GEANT4_USE_FREETYPE)
  find_package(Freetype REQUIRED)
endif()

GEANT4_ADD_FEATURE(GEANT4_USE_FREETYPE "Building Geant4 analysis library with Freetype support")

