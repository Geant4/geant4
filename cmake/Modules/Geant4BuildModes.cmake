# - Add Geant4 specific build modes and additional flags
#
# This follows the guide on adding a new mode on the CMake wiki:
#
# http://www.cmake.org/Wiki/CMake_FAQ#How_can_I_extend_the_build_modes_with_a_custom_made_one_.3F
#
# Geant4 supports the standard CMake build types/configurations of
#
# Release Debug MinSizeRel RelWithDebInfo
#
# In addition, two types specifically for development are added:
#
# TestRelease:
#   For trial production and extended testing. It has verbose
#   output, has debugging symbols, and adds definitions to allow FPE
#   and physics conservation law testing where supported.
#
# Maintainer:
#   For development of the toolkit. It adds debugging, and enables the use
#   of library specific debugging via standardized definitions.
#
# Compiler flags specific to these build types are set in the cache, and the 
# types are added to the CMAKE_BUILD_TYPE cache string and to
# CMAKE_CONFIGURATION_TYPES if appropriate to the build tool being used.
#

#------------------------------------------------------------------------------
# Add TestRelease{Debug} Modes and cache init flags
#
set(CMAKE_CXX_FLAGS_TESTRELEASE "${CMAKE_CXX_FLAGS_TESTRELEASE_INIT}"
  CACHE STRING "Flags used by the compiler during TestRelease builds"
)


#------------------------------------------------------------------------------
# Add Maintainer Mode
#
set(CMAKE_CXX_FLAGS_MAINTAINER "${CMAKE_CXX_FLAGS_MAINTAINER_INIT}"
  CACHE STRING "Flags used by the compiler during Maintainer builds"
)

#----------------------------------------------------------------------------
# Mark all the additional mode flags as advanced because most users will 
# never need to see them
mark_as_advanced(
  CMAKE_CXX_FLAGS_TESTRELEASE
  CMAKE_CXX_FLAGS_MAINTAINER
)

#------------------------------------------------------------------------------
# Add the new configuration types ONLY if the build tool supports multiple
# configurations
#
if(CMAKE_CONFIGURATION_TYPES)
  list(APPEND CMAKE_CONFIGURATION_TYPES TestRelease)
  list(APPEND CMAKE_CONFIGURATION_TYPES Maintainer)
  list(REMOVE_DUPLICATES CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" 
    CACHE STRING "Geant4 configurations for multimode build tools"
    FORCE
    )
endif()

#------------------------------------------------------------------------------
# Update build type information ONLY for single mode build tools, adding
# default type if none has been set, but otherwise leaving value alone.
# NB: this doesn't allow "None" for the build type - would need something
# more sophiticated using an internal cache variable.
# Good enough for now!
#
if(NOT CMAKE_CONFIGURATION_TYPES)
  if(NOT CMAKE_BUILD_TYPE)
    # Default to a Release build if nothing else...
    set(CMAKE_BUILD_TYPE Release
      CACHE STRING "Choose the type of build, options are: None Release TestRelease MinSizeRel Debug RelWithDebInfo MinSizeRel Maintainer."
      FORCE
      )
  else()
    # Force to the cache, but use existing value.
    set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}"
      CACHE STRING "Choose the type of build, options are: None Release TestRelease MinSizeRel Debug RelWithDebInfo MinSizeRel Maintainer."
      FORCE
      )
  endif()
endif()
