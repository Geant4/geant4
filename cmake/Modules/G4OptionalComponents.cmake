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
#  USOLIDS  - Allow use of VecGeom classes in geometry, and must find
#             external VecGeom install
#  FREETYPE - For analysis/vis modules, find external FreeType install

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
set(_default_use_system_clhep OFF)
if(CLHEP_ROOT_DIR)
  set(_default_use_system_clhep ON)
  list(INSERT CMAKE_PREFIX_PATH 0 "${CLHEP_ROOT_DIR}")
endif()

option(GEANT4_USE_SYSTEM_CLHEP "Use system CLHEP library" ${_default_use_system_clhep})

if(GEANT4_USE_SYSTEM_CLHEP)
  find_package(CLHEP 2.4.6.0 REQUIRED CONFIG)
  geant4_save_package_variables(CLHEP CLHEP_DIR)
else()
  set(CLHEP_FOUND TRUE)
  set(CLHEP_LIBRARIES G4clhep)
endif()

geant4_add_feature(GEANT4_USE_SYSTEM_CLHEP "Using system CLHEP library")

#-----------------------------------------------------------------------
# Find required EXPAT package
# We always use the internal G4expat on WIN32.
# On other platforms, we default to use the system library.
# If we use the internal G4expat, fix the EXPAT_XXX variables to point
# to the internal location of expat headers and library target so Geant4
# clients of expat have a consistent interface.
#
if(WIN32)
  set(EXPAT_FOUND TRUE)
  set(GEANT4_USE_BUILTIN_EXPAT TRUE)
  set(EXPAT_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/source/externals/expat/include)
  set(EXPAT_LIBRARIES G4expat)
else()
  option(GEANT4_USE_SYSTEM_EXPAT "Use system Expat library" ON)

  if(GEANT4_USE_SYSTEM_EXPAT)
    # If system package requested, find it or fail
    find_package(EXPAT REQUIRED)
    # Shim only needed until we require CMake >= 3.10
    include("${CMAKE_CURRENT_LIST_DIR}/G4EXPATShim.cmake")

    # Check version requirement externally to provide information
    # on using internal expat.
    if(${EXPAT_VERSION_STRING} VERSION_LESS "2.0.1")
      set(__badexpat_include_dir ${EXPAT_INCLUDE_DIR})
      set(__badexpat_library ${EXPAT_LIBRARY})
      unset(EXPAT_FOUND)
      unset(EXPAT_INCLUDE_DIR CACHE)
      unset(EXPAT_LIBRARY CACHE)
      message(FATAL_ERROR
"Detected system expat header and library:
EXPAT_INCLUDE_DIR = ${__badexpat_include_dir}
EXPAT_LIBRARY = ${__badexpat_library}
are of insufficient version '${EXPAT_VERSION_STRING}' (Required >= 2.0.1)
Set the above CMake variables to point to an expat install of the required version, or set GEANT4_USE_SYSTEM_EXPAT to OFF to use Geant4's packaged version.")
    endif()

    # Backward compatibility for sources.cmake using the variable
    set(EXPAT_LIBRARIES EXPAT::EXPAT)
    geant4_save_package_variables(EXPAT EXPAT_INCLUDE_DIR EXPAT_LIBRARY)
  else()
    set(EXPAT_FOUND TRUE)
    set(GEANT4_USE_BUILTIN_EXPAT TRUE)
    set(EXPAT_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/source/externals/expat/include)
    set(EXPAT_LIBRARIES G4expat)
  endif()
endif()

geant4_add_feature(GEANT4_USE_SYSTEM_EXPAT "Using system EXPAT library")

#-----------------------------------------------------------------------
# Find required ZLIB package, defaulting in internal
option(GEANT4_USE_SYSTEM_ZLIB "Use system zlib library" OFF)
if(GEANT4_USE_SYSTEM_ZLIB)
  find_package(ZLIB REQUIRED)
  # Backward compatibility for sources.cmake using the variable
  set(ZLIB_LIBRARIES ZLIB::ZLIB)
  geant4_save_package_variables(ZLIB ZLIB_INCLUDE_DIR ZLIB_LIBRARY_DEBUG ZLIB_LIBRARY_RELEASE)
else()
  set(ZLIB_FOUND TRUE)
  set(GEANT4_USE_BUILTIN_ZLIB TRUE)
  set(ZLIB_LIBRARIES G4zlib)
endif()

geant4_add_feature(GEANT4_USE_SYSTEM_ZLIB "Using system zlib library")

#-----------------------------------------------------------------------
# Find required PTL, optional TBB support, defaulting to internal
option(GEANT4_USE_SYSTEM_PTL "Use system PTL library" OFF)
option(GEANT4_USE_TBB "Use TBB as PTL's tasking backend" OFF)
cmake_dependent_option(GEANT4_USE_PTL_LOCKS "Enable mutex locking in PTL task subqueues" OFF "NOT GEANT4_USE_SYSTEM_PTL" OFF)
mark_as_advanced(GEANT4_USE_SYSTEM_PTL GEANT4_USE_PTL_LOCKS GEANT4_USE_TBB)

if(GEANT4_USE_SYSTEM_PTL)
  if(GEANT4_USE_TBB)
    find_package(PTL 2.0.0 REQUIRED TBB)
  else()
    find_package(PTL 2.0.0 REQUIRED)
    if(PTL_TBB_FOUND AND NOT GEANT4_USE_TBB)
      set(GEANT4_USE_TBB ON CACHE BOOL "Use TBB as PTL's tasking backend (forced)" FORCE)
      message(STATUS "Forcing GEANT4_USE_TBB to ON, required by system PTL: ${PTL_DIR}")
    endif()
  endif() 
  set(PTL_LIBRARIES PTL::ptl)
  geant4_save_package_variables(PTL PTL_DIR)
else()
  # Prempt PTL's find of TBB so we can cache the variables
  if(GEANT4_USE_TBB)
    find_package(TBB REQUIRED)
  endif()
  set(PTL_FOUND TRUE)
  set(GEANT4_USE_BUILTIN_PTL TRUE)
  set(PTL_LIBRARIES G4ptl)
endif()

geant4_add_feature(GEANT4_USE_SYSTEM_PTL "Using system PTL library")
geant4_add_feature(GEANT4_USE_PTL_LOCKS "Enabled mutex locking in PTL task subqueues")
geant4_add_feature(GEANT4_USE_TBB "Using TBB as PTL's tasking backend")
# Always cache TBB even if not used. Empty vars will not be exported.
geant4_save_package_variables(TBB TBB_DIR TBB_INCLUDE_DIR TBB_LIBRARY TBB_LIBRARY_DEBUG TBB_LIBRARY_RELEASE)

#-----------------------------------------------------------------------
# Optional Support for GDML - requires Xerces-C package
# Relies on XercesC::XercesC imported target (since CMake 3.5)
set(_default_use_gdml OFF)
if(XERCESC_ROOT_DIR)
  set(_default_use_gdml ON)
  list(INSERT CMAKE_PREFIX_PATH 0 "${XERCESC_ROOT_DIR}")
endif()

option(GEANT4_USE_GDML "Build Geant4 with GDML support" ${_default_use_gdml})

if(GEANT4_USE_GDML)
  find_package(XercesC REQUIRED)
  geant4_save_package_variables(XercesC XercesC_INCLUDE_DIR XercesC_LIBRARY_DEBUG XercesC_LIBRARY_RELEASE)
endif()

geant4_add_feature(GEANT4_USE_GDML "Building Geant4 with GDML support")

#-----------------------------------------------------------------------
# Optional use of smart stack
#  With this option, G4StackManager uses G4SmartTrackStack instead of
# ordinary G4TrackStack as the Urgent stack. G4SmartTrackStack tries to
# stick to the same kind of particle as the previous track when Pop()
# is called. This G4SmartTrackStack may provide some performance
# improvements in particular for crystal calorimeters in high energy
# physics experiments. On the other hand, G4SmartTrackStack won't give
# any benefit for granular geometry or lower energy applications, while
# it may causes some visible memory footprint increase.

option(GEANT4_USE_SMARTSTACK "Use smart track stack" OFF)
mark_as_advanced(GEANT4_USE_SMARTSTACK)
geant4_add_feature(GEANT4_USE_SMARTSTACK "Use smart track stack")

#-----------------------------------------------------------------------
# Optional Support for TiMemory -- timing, memory, HW counters, roofline, gperftools, etc.
# easily installed via:
#   git clone https://github.com/NERSC/timemory.git timemory
#   pip install -vvv ./timemory
# and prepending CMAKE_PREFIX_PATH with `python -c "import sys; print(sys.prefix)"`
#
set(_default_use_timemory OFF)
if(TiMemory_DIR)
  set(_default_use_timemory ON)
endif()

option(GEANT4_USE_TIMEMORY "Build Geant4 with TiMemory support" ${_default_use_timemory})
mark_as_advanced(GEANT4_USE_TIMEMORY)

if(GEANT4_USE_TIMEMORY)
  # by default just use the library which will import all the components that it
  # was built with but once can add more
  if(BUILD_SHARED_LIBS)
    set(_G4timemory_DEFAULT_COMPONENTS cxx shared OPTIONAL_COMPONENTS)
  else()
    set(_G4timemory_DEFAULT_COMPONENTS cxx static OPTIONAL_COMPONENTS)
  endif()
  set(G4timemory_COMPONENTS "${_G4timemory_DEFAULT_COMPONENTS}" CACHE STRING
      "timemory INTERFACE libraries that activate various capabilities in toolkit")
  set(G4timemory_VERSION 3.2)
  set(timemory_FIND_COMPONENTS_INTERFACE geant4-timemory)
  add_library(geant4-timemory INTERFACE)
  find_package(timemory ${G4timemory_VERSION} REQUIRED COMPONENTS ${G4timemory_COMPONENTS})
  install(TARGETS geant4-timemory EXPORT Geant4LibraryDepends)
  set(timemory_LIBRARIES geant4-timemory)
  geant4_save_package_variables(timemory timemory_DIR)
endif()

geant4_add_feature(GEANT4_USE_TIMEMORY "Building Geant4 with TiMemory support")

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
# Optional support for VecGeom - Requires VecGeom package
# GEANT4_USE_USOLIDS is a list argument which must take values:
# - Empty/CMake boolean false : DEFAULT, do not replace G4 solids with VecGeom
# - ALL/CMake boolean true : Replace all G4 solids with VecGeom
# - List : Replace listed G4 solids with VecGeom equivalents. Valid elements are
set(GEANT4_USOLIDS_SHAPES
  BOX
  CONS
  CTUBS
  ELLIPSOID
  ELLIPTICALCONE
  ELLIPTICALTUBE
  EXTRUDEDSOLID
  HYPE
  GENERICPOLYCONE
  GENERICTRAP
  ORB
  PARA
  PARABOLOID
  POLYCONE
  POLYHEDRA
  SPHERE
  TESSELLATEDSOLID
  TET
  TRAP
  TRD
  TORUS
  TUBS
  )

set(GEANT4_USE_USOLIDS OFF CACHE STRING "EXPERIMENTAL: List Geant4 solids to replace with VecGeom equivalents (ALL;${GEANT4_USOLIDS_SHAPES})")
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


# - Geant4 USolids/VecGom setup
if(GEANT4_USE_ALL_USOLIDS OR GEANT4_USE_PARTIAL_USOLIDS)
  find_package(VecGeom 1.2.0 REQUIRED)
  geant4_save_package_variables(VecGeom VecGeom_DIR)

  # If VecCore_DIR is set, means updated VecGeom install used, so
  # also store VecCore_DIR
  if(VecCore_DIR)
    geant4_save_package_variables(VecGeom VecCore_DIR)
  endif()

  if(GEANT4_USE_ALL_USOLIDS)
    set(G4GEOM_USE_USOLIDS TRUE)
    geant4_add_feature(GEANT4_USE_USOLIDS "Replacing Geant4 solids with all VecGeom equivalents (EXPERIMENTAL)")
  else()
    set(G4GEOM_USE_PARTIAL_USOLIDS TRUE)
    foreach(__g4_usolid_shape ${GEANT4_USE_PARTIAL_USOLIDS_SHAPE_LIST})
      set(G4GEOM_USE_U${__g4_usolid_shape} TRUE)
    endforeach()
    geant4_add_feature(GEANT4_USE_USOLIDS "Replacing Geant4 solids with VecGeom equivalents for ${GEANT4_USE_PARTIAL_USOLIDS_SHAPE_LIST} (EXPERIMENTAL)")
  endif()
endif()


#-----------------------------------------------------------------------
# Optional support for Freetype - Requires external Freetype install
#
option(GEANT4_USE_FREETYPE "Build Geant4 Analysis/Visualization with Freetype support" OFF)
mark_as_advanced(GEANT4_USE_FREETYPE)

if(GEANT4_USE_FREETYPE)
  find_package(Freetype REQUIRED)
  geant4_save_package_variables(Freetype
    FREETYPE_INCLUDE_DIR_freetype2
    FREETYPE_INCLUDE_DIR_ft2build
    FREETYPE_LIBRARY_DEBUG
    FREETYPE_LIBRARY_RELEASE)
endif()

geant4_add_feature(GEANT4_USE_FREETYPE "Building Geant4 Analysis/Visualization with Freetype support")

#-----------------------------------------------------------------------
# Optional support for HDF5
# - Requires external HDF5 1.8 or higher install
# - Install must be MT safe if building Geant4 in MT mode
#
option(GEANT4_USE_HDF5 "Build Geant4 analysis library with HDF5 support" OFF)
mark_as_advanced(GEANT4_USE_HDF5)

if(GEANT4_USE_HDF5)
  find_package(HDF5 1.8 REQUIRED)
  include("${CMAKE_CURRENT_LIST_DIR}/G4HDF5Shim.cmake")
  # Backward compatibility
  set(HDF5_LIBRARIES hdf5::hdf5)

  # May have found via config mode...
  if(HDF5_DIR)
    geant4_save_package_variables(HDF5 HDF5_DIR)
  else()
    # Otherwise almost certainly used compiler wrapper
    geant4_save_package_variables(HDF5 HDF5_C_COMPILER_EXECUTABLE HDF5_C_LIBRARY_hdf5)
  endif()
endif()

geant4_add_feature(GEANT4_USE_HDF5 "Building Geant4 analysis library with HDF5 support")
