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
# Find required CLHEP package
# We prefer to use our internal CLHEP library, but provide an option
# to allow an external CLHEP to be used should the user require this.
# We also allow that it can be automatically enabled by providing
# the CLHEP_ROOT_DIR option (which FindCLHEP will recognize)
#
# As requested by ATLAS, an additional option for preferring use of
# CLHEP's granular libs is provided when using a system CLHEP.
#
<<<<<<< HEAD
# KNOWNISSUE : For internal CLHEP, how to deal with static and shared?
if(CLHEP_ROOT_DIR)
  set(_default_use_system_clhep ON)
else()
  set(_default_use_system_clhep OFF)
=======
set(_default_use_system_clhep OFF)
if(CLHEP_ROOT_DIR)
  set(_default_use_system_clhep ON)
  list(INSERT CMAKE_PREFIX_PATH 0 "${CLHEP_ROOT_DIR}")
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
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

<<<<<<< HEAD
  # Find CLHEP using package-mode (i.e. FindCLHEP)
  # This will set up imported targets for us, but we need
  # to set CLHEP_LIBRARIES afterwards to use these
  find_package(CLHEP 2.3.1.0 REQUIRED ${__g4_clhep_components})

  set(CLHEP_LIBRARIES CLHEP::CLHEP)
  if(GEANT4_USE_SYSTEM_CLHEP_GRANULAR)
    set(CLHEP_LIBRARIES)
    foreach(_comp ${__g4_clhep_components})
      list(APPEND CLHEP_LIBRARIES  "CLHEP::${_comp}")
    endforeach()
  endif()

  set(GEANT4_USE_SYSTEM_CLHEP TRUE)
=======
  find_package(CLHEP 2.3.3.0 REQUIRED ${__g4_clhep_components} CONFIG)

  geant4_save_package_variables(CLHEP CLHEP_DIR)
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
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
<<<<<<< HEAD
=======
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
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
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
# Rely on ZLIB::ZLIB imported target (since CMake 3.1)
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
<<<<<<< HEAD
=======
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
#       and setting timemory_DIR to `python -c "import sys; print(sys.prefix)"`
#
set(_default_use_timemory OFF)
if(TiMemory_DIR)
  set(_default_use_timemory ON)
endif()

option(GEANT4_USE_TIMEMORY "Build Geant4 with TiMemory support" ${_default_use_timemory})
mark_as_advanced(GEANT4_USE_TIMEMORY)

if(GEANT4_USE_TIMEMORY)
  set(_G4timemory_DEFAULT_COMPONENTS headers caliper papi gotcha gperftools-cpu vector)
  set(G4timemory_COMPONENTS "${_G4timemory_DEFAULT_COMPONENTS}" CACHE STRING
      "timemory INTERFACE libraries that activate various capabilities in toolkit")
  set(G4timemory_VERSION 3.0)
  set(timemory_FIND_COMPONENTS_INTERFACE geant4-timemory)
  find_package(timemory ${G4timemory_VERSION} REQUIRED COMPONENTS ${G4timemory_COMPONENTS})
  set(timemory_LIBRARIES geant4-timemory)
  geant4_save_package_variables(timemory timemory_DIR)
endif()

geant4_add_feature(GEANT4_USE_TIMEMORY "Building Geant4 with TiMemory support")

#-----------------------------------------------------------------------
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
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
<<<<<<< HEAD
=======
  CTUBS
  ELLIPSOID
  ELLIPTICALCONE
  ELLIPTICALTUBE
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
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
<<<<<<< HEAD
  find_package(USolids REQUIRED)

  if(GEANT4_USE_ALL_USOLIDS)
    set(GEANT4_USOLIDS_COMPILE_DEFINITIONS "-DG4GEOM_USE_USOLIDS")
    GEANT4_ADD_FEATURE(GEANT4_USE_USOLIDS "Replacing all Geant4 solids with USolids equivalents (EXPERIMENTAL)")
=======
  # VecGeom's config file doesn't support versioning...
  find_package(VecGeom REQUIRED)
  # Shim until VecGeom supports config mode properly
  include("${CMAKE_CURRENT_LIST_DIR}/G4VecGeomShim.cmake")
  # Backward Compatibility
  set(VECGEOM_LIBRARIES VecGeom::VecGeom)

  geant4_save_package_variables(VecGeom VecGeom_DIR)

  if(GEANT4_USE_ALL_USOLIDS)
    set(G4GEOM_USE_USOLIDS TRUE)
    GEANT4_ADD_FEATURE(GEANT4_USE_USOLIDS "Replacing Geant4 solids with all VecGeom equivalents (EXPERIMENTAL)")
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
  else()
    set(G4GEOM_USE_PARTIAL_USOLIDS TRUE)
    foreach(__g4_usolid_shape ${GEANT4_USE_PARTIAL_USOLIDS_SHAPE_LIST})
      set(G4GEOM_USE_U${__g4_usolid_shape} TRUE)
    endforeach()
    GEANT4_ADD_FEATURE(GEANT4_USE_USOLIDS "Replacing Geant4 solids with USolids equivalents for ${GEANT4_USE_PARTIAL_USOLIDS_SHAPE_LIST} (EXPERIMENTAL)")
  endif()
<<<<<<< HEAD

  # Combined definitions
  add_definitions(${GEANT4_USOLIDS_COMPILE_DEFINITIONS})

  # Add USolids inc dirs here - can be removed once USolids supports
  # INTERFACE_INCLUDE_DIRECTORIES
  include_directories(${USOLIDS_INCLUDE_DIRS})
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
endif()


#-----------------------------------------------------------------------
# Optional support for Freetype - Requires external Freetype install
#
option(GEANT4_USE_FREETYPE "Build Geant4 analysis library with Freetype support" OFF)
mark_as_advanced(GEANT4_USE_FREETYPE)

if(GEANT4_USE_FREETYPE)
  find_package(Freetype REQUIRED)
  # Shim only needed until we require CMake >= 3.10
  include("${CMAKE_CURRENT_LIST_DIR}/G4FreetypeShim.cmake")

  geant4_save_package_variables(Freetype
    FREETYPE_INCLUDE_DIR_freetype2
    FREETYPE_INCLUDE_DIR_ft2build
    FREETYPE_LIBRARY_DEBUG
    FREETYPE_LIBRARY_RELEASE)
endif()

geant4_add_feature(GEANT4_USE_FREETYPE "Building Geant4 analysis library with Freetype support")

<<<<<<< HEAD
=======
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
  set(HDF5_LIBRARIES Geant4::HDF5)

  # May have found via config mode...
  if(HDF5_DIR)
    geant4_save_package_variables(HDF5 HDF5_DIR)
  else()
    # Otherwise almost certainly used compiler wrapper
    geant4_save_package_variables(HDF5
      HDF5_C_COMPILER_EXECUTABLE
      HDF5_C_LIBRARY_hdf5)
    endif()
endif()

GEANT4_ADD_FEATURE(GEANT4_USE_HDF5 "Building Geant4 analysis library with HDF5 support")

>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
