# - Use file for Geant4
# This file should be included after a find_package call has successfully
# located Geant4. If Geant4 has been located via the Geant4Config.cmake
# config file, this will have set the following variable:
#
#  Geant4_USE_FILE : Point to the location of the use file for the found
#                    Geant4 installation.
#
# Inclusion of this file, e.g. via
#
#  include(${Geant4_USE_FILE})
#
# results in the addition of the Geant4 compile definitions and
# include directories to those of the directory in which this file is
# included.
#
# Header paths are added to include_directories as SYSTEM type directories
#
# The recommended Geant4 compiler flags are also prepended to
# CMAKE_CXX_FLAGS but duplicated flags are NOT removed. This permits
# client of UseGeant4 to override Geant4's recommended flags if required
# and at their own risk.
#
# Advanced users requiring special sets of flags, or the removal of
# duplicate flags should therefore *not* use this file, preferring the
# direct use of the Geant4_XXXX variables set by the Geant4Config file.
#
# The last thing the module does is to optionally include an internal Use
# file. This file can contain variables, functions and macros for strict
# internal use in Geant4, such as building and running validation tests.
#

#-----------------------------------------------------------------------
# We need to set the compile definitions and include directories
#
add_definitions(${Geant4_DEFINITIONS})
include_directories(AFTER SYSTEM ${Geant4_INCLUDE_DIRS})

#-----------------------------------------------------------------------
# Because Geant4 is sensitive to the compiler flags, let's set the base
# set here. This reproduces as far as possible the behaviour of the
# original makefile system. However, we append any existing CMake flags in
# case the user wishes to override these (at their own risk).
# Though this may lead to duplication, that should not affect behaviour.
#
set(CMAKE_CXX_FLAGS                "${Geant4_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "${Geant4_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${Geant4_CXX_FLAGS_MINSIZEREL} ${CMAKE_CXX_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE        "${Geant4_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${Geant4_CXX_FLAGS_RELWITHDEBINFO} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set(CMAKE_EXE_LINKER_FLAGS         "${Geant4_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")

#-----------------------------------------------------------------------
# Locate ourselves
#
get_filename_component(_use_geant4_dir ${CMAKE_CURRENT_LIST_FILE} PATH)

#-----------------------------------------------------------------------
# Append the local module path to CMAKE_MODULE_PATH to automatically
# make FindXXX modules for examples available
#
list(APPEND CMAKE_MODULE_PATH ${_use_geant4_dir}/Modules)

#-----------------------------------------------------------------------
# Include internal use file if it exists. It should only exist in the
# build tree!
#
include(${_use_geant4_dir}/UseGeant4_internal.cmake OPTIONAL)

