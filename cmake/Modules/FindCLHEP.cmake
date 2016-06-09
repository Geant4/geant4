# - Try to find the CLHEP High Energy Physics library and headers
# Usage of this module is as follows
#
# == Using any header-only components of CLHEP: ==
#
#     find_package( CLHEP 2.0.4.2 )
#     if(CLHEP_FOUND)
#         include_directories(${CLHEP_INCLUDE_DIRS})
#         add_executable(foo foo.cc)
#     endif()
#
# == Using the binary CLHEP library ==
#
#     find_package( CLHEP 2.0.4.2 )
#     if(CLHEP_FOUND)
#         include_directories(${CLHEP_INCLUDE_DIRS})
#         add_executable(foo foo.cc)
#         target_link_libraries(foo ${CLHEP_LIBRARIES})
#     endif()
#
# You can provide a minimum version number that should be used.
# If you provide this version number and specify the REQUIRED attribute,
# this module will fail if it can't find a CLHEP of the specified version
# or higher. If you further specify the EXACT attribute, then this module 
# will fail if it can't find a CLHEP with a version eaxctly as specified.
#
# ===========================================================================
# Variables used by this module which can be used to change the default
# behaviour, and hence need to be set before calling find_package:
#
#  CLHEP_ROOT_DIR        The preferred installation prefix for searching for 
#                        CLHEP. Set this if the module has problems finding 
#                        the proper CLHEP installation.
#
# If you don't supply CLHEP_ROOT_DIR, the module will search on the standard
# system paths. On UNIX, the module will also try to find the clhep-config
# program in the PATH, and if found will use the prefix supplied by this
# program as a HINT on where to find the CLHEP headers and libraries.
#
# You can re-run CMake with a different version of CLHEP_ROOT_DIR to 
# force a new search for CLHEP using the new version of CLHEP_ROOT_DIR.
#
# ============================================================================
# Variables set by this module:
#
#  CLHEP_FOUND           System has CLHEP.
#
#  CLHEP_INCLUDE_DIRS    CLHEP include directories: not cached.
#
#  CLHEP_LIBRARIES       Link to these to use the CLHEP library: not cached.
#
# ===========================================================================
# If CLHEP is installed in a non-standard way, e.g. a non GNU-style install
# of <prefix>/{lib,include}, then this module may fail to locate the headers
# and libraries as needed. In this case, the following cached variables can
# be editted to point to the correct locations.
#
#  CLHEP_INCLUDE_DIR    The path to the CLHEP include directory: cached
#
#  CLHEP_LIBRARY        The path to the CLHEP library: cached
# 
# You should not need to set these in the vast majority of cases
#

#============================================================================
# Copyright (C) 2010,2011 Ben Morgan <Ben.Morgan@warwick.ac.uk>
# Copyright (C) 2010,2011 University of Warwick
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, 
#   this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice, 
#   this list of conditions and the following disclaimer in the documentation 
#   and/or other materials provided with the distribution.
#
# * Neither the name of the University of Warwick nor the names of its 
#   contributors may be used to endorse or promote products derived from 
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#============================================================================


#----------------------------------------------------------------------------
# Enable re-search if known CLHEP_ROOT_DIR changes?
#
if(NOT "${CLHEP_ROOT_DIR}" STREQUAL "${CLHEP_INTERNAL_ROOT_DIR}")
    if(CLHEP_INTERNAL_ROOT_DIR AND NOT CLHEP_FIND_QUIETLY)
        message(STATUS "CLHEP_ROOT_DIR Changed, Rechecking for CLHEP")
    endif()

    set(CLHEP_INTERNAL_ROOT_DIR ${CLHEP_ROOT_DIR} 
        CACHE INTERNAL "Last value supplied for where to locate CLHEP")
    set(CLHEP_INCLUDE_DIR CLHEP_INCLUDE_DIR-NOTFOUND)
    set(CLHEP_LIBRARY CLHEP_LIBRARY-NOTFOUND)
    set(CLHEP_CONFIG_EXECUTABLE CLHEP_CONFIG_EXECUTABLE-NOTFOUND)
    set(CLHEP_LIBRARIES )
    set(CLHEP_INCLUDE_DIRS )
    set(CLHEP_FOUND FALSE)
endif()

#----------------------------------------------------------------------------
# - If we already found CLHEP, be quiet
#
if(CLHEP_INCLUDE_DIR AND CLHEP_LIBRARY)
    set(CLHEP_FIND_QUIETLY TRUE)
endif()

#----------------------------------------------------------------------------
# Set up HINTS on where to look for CLHEP
# If we're on UNIX, see if we can find clhep-config and use its --prefix
# as an extra hint.
#
set(_clhep_root_hints ${CLHEP_ROOT_DIR})

if(UNIX)
    # Try and find clhep-config in the user's path, but hint at the bin
    # directory under CLHEP_ROOT_DIR because we'd ideally like to pick up
    # the config program that matches the libraries/headers.
    # We only use it as a fallback though.
    find_program(CLHEP_CONFIG_EXECUTABLE clhep-config
        HINTS ${_clhep_root_hints}/bin
        DOC "Path to CLHEP's clhep-config program")
    mark_as_advanced(CLHEP_CONFIG_EXECUTABLE)

    if(CLHEP_CONFIG_EXECUTABLE)
        execute_process(COMMAND ${CLHEP_CONFIG_EXECUTABLE} --prefix
            OUTPUT_VARIABLE _clhep_config_prefix
            OUTPUT_STRIP_TRAILING_WHITESPACE)

        list(APPEND _clhep_root_hints ${_clhep_config_prefix})
    endif()
elseif(WIN32 AND NOT UNIX)
    # Do we need to set suitable defaults?
endif()


#----------------------------------------------------------------------------
# Find the CLHEP headers
# Use Units/defs.h as locator as this is pretty consistent through versions
find_path(CLHEP_INCLUDE_DIR CLHEP/Units/defs.h
    HINTS ${_clhep_root_hints}
    PATH_SUFFIXES include
    DOC "Path to the CLHEP headers" 
)

#----------------------------------------------------------------------------
# Find the CLHEP library
# Prefer lib64 if available.
find_library(CLHEP_LIBRARY CLHEP
    HINTS ${_clhep_root_hints}
    PATH_SUFFIXES lib64 lib
    DOC "Path to the CLHEP library"
)


#----------------------------------------------------------------------------
# Extract the CLHEP version from defs.h
# Versions COMPATIBLE if RequestedVersion > FoundVersion
# Also check if versions exact

if(CLHEP_INCLUDE_DIR)
    set(CLHEP_VERSION 0)
    file(READ "${CLHEP_INCLUDE_DIR}/CLHEP/Units/defs.h" _CLHEP_DEFS_CONTENTS)
    string(REGEX REPLACE ".*#define PACKAGE_VERSION \"([0-9.]+).*" "\\1"
        CLHEP_VERSION "${_CLHEP_DEFS_CONTENTS}")

    if(NOT CLHEP_FIND_QUIETLY)
        message(STATUS "Found CLHEP Version ${CLHEP_VERSION}")
    endif()
  
    if(CLHEP_FIND_VERSION)
        set(CLHEP_VERSIONING_TESTS CLHEP_VERSION_COMPATIBLE)

        if("${CLHEP_VERSION}" VERSION_LESS "${CLHEP_FIND_VERSION}")
            set(CLHEP_VERSION_COMPATIBLE FALSE)
        else()
            set(CLHEP_VERSION_COMPATIBLE TRUE)

            if(CLHEP_FIND_VERSION_EXACT)
                if("${CLHEP_VERSION}" VERSION_EQUAL "${CLHEP_FIND_VERSION}")
                    set(CLHEP_VERSION_EXACT TRUE)
                endif()
                    list(APPEND CLHEP_VERSIONING_TESTS CLHEP_VERSION_EXACT)
            endif()
        endif()
    endif()
endif()

#----------------------------------------------------------------------------
# Construct an error message for FPHSA
#
set(CLHEP_DEFAULT_MSG "Could NOT find CLHEP:\n")

if(NOT CLHEP_INCLUDE_DIR)
    set(CLHEP_DEFAULT_MSG "${CLHEP_DEFAULT_MSG}CLHEP Header Path Not Found\n")
endif()

if(NOT CLHEP_LIBRARY)
    set(CLHEP_DEFAULT_MSG "${CLHEP_DEFAULT_MSG}CLHEP Library Not Found\n")
endif()

if(CLHEP_FIND_VERSION)
    if(NOT CLHEP_VERSION_COMPATIBLE)
        set(CLHEP_DEFAULT_MSG "${CLHEP_DEFAULT_MSG}Incompatible versions, ${CLHEP_VERSION}(found) < ${CLHEP_FIND_VERSION}(required)\n")
    endif()

    if(CLHEP_FIND_VERSION_EXACT)
        if(NOT CLHEP_VERSION_EXACT)
            set(CLHEP_DEFAULT_MSG "${CLHEP_DEFAULT_MSG}Non-exact versions, ${CLHEP_VERSION}(found) != ${CLHEP_FIND_VERSION}(required)\n")
        endif()
    endif()
endif()




#----------------------------------------------------------------------------
# Handle the QUIETLY and REQUIRED arguments, setting CLHEP_FOUND to TRUE if
# all listed variables are TRUE
#
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CLHEP 
    "${CLHEP_DEFAULT_MSG}"
    CLHEP_LIBRARY
    CLHEP_INCLUDE_DIR
    ${CLHEP_VERSIONING_TESTS}
    )


#----------------------------------------------------------------------------
# If we found CLHEP, set the needed non-cache variables
#
if(CLHEP_FOUND)
    set(CLHEP_LIBRARIES ${CLHEP_LIBRARY})
    set(CLHEP_INCLUDE_DIRS ${CLHEP_INCLUDE_DIR})
endif()

#----------------------------------------------------------------------------
# Mark cache variables that can be adjusted as advanced
#
mark_as_advanced(CLHEP_INCLUDE_DIR CLHEP_LIBRARY)

