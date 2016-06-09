# - Try to find CLHEP library on *NIX Platforms
#
# Once done this will define
#
#   CLHEP_FOUND - system has the required CLHEP library
#   CLHEP_INCLUDE_DIRS- the CLHEP include directory
#   CLHEP_LIBRARIES - the libraries needed to use CLHEP
#   CLHEP_VERSION_OK - Flag indicating if discovered version satisfies
#                      versioning requirements.
#
# Usage of this module as follows:
#
#   FIND_PACKAGE(CLHEP 2.0.3.3 [EXACT] [REQUIRED|QUIET])
#
# Variables used by this module, which can change the default behaviour and
# need to be set before calling find_package:
#
#
#   CLHEP_CONFIG_EXECUTABLE     All information on the CLHEP installation is
#                               derived from the clhep-config program. Set this
#                               variable to the full path to clhep-config if 
#                               the module has problems finding the proper 
#                               CLHEP installation (Cached).
# 
#   CLHEP_MAX_VERSION           Set this variable to maximum version number
#                               required, in the format "x.y.z" or "x.y". This
#                               is used with any version number supplied to
#                               FIND_PACKAGE to do a range check on the
#                               discovered version.
#
# The module uses the following logic to perform version checking
#
# 1) FIND_PACKAGE version a.b.c supplied, EXACT set:
#   Check and require that discovered version x.y.z == a.b.c
#   CLHEP_MAX_VERSION ignored.
#
# 2) FIND_PACKAGE version a.b.c supplied
#   i) If CLHEP_MAX_VERSION not set, check and require that discovered version
#      x.y.z >= a.b.c
#   ii) If CLHEP_MAX_VERSION p.q.r set, check and require that discovered
#       version x.y.z satisfies a.b.c <= x.y.z < p.q.r
#
# 3) FIND_PACKAGE version a.b.c NOT supplied
#   i) If CLHEP_MAX_VERSION not set, not version check performed.
#   ii) If CLHEP_MAX_VERSION p.q.r set, check and require that discovered
#       version x.y.z < p.q.r
#

# Copyright (c) 2009, Ben Morgan, <Ben.Morgan@warwick.ac.uk>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.


include(ResolveCompilerPaths)
include(MacroEnsureVersion)

if (CLHEP_FOUND)
    # We already found CLHEP
    # For now ignore different versioning between calls - yes makes no sense
    # for subcalls to recheck?
    set(CLHEP_FIND_QUIETLY TRUE)

else ()
    # Find the clhep-config program
    find_program(CLHEP_CONFIG_EXECUTABLE NAMES clhep-config)

    # Reset vars
    set(CLHEP_LIBRARIES)
    set(CLHEP_INCLUDE_DIRS)
    set(CLHEP_VERSION_OK)

    # If we find it...
    if (CLHEP_CONFIG_EXECUTABLE)
        # Check version first
        execute_process(COMMAND ${CLHEP_CONFIG_EXECUTABLE} --version OUTPUT_VARIABLE CLHEP_DISCOVERED_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "Found CLHEP version: ${CLHEP_DISCOVERED_VERSION}")

        # VERSIONING CHECK
        # At the end of this section, CLHEP_VERSION_OK will be true if
        # discovered version satisfies requirements
        if (CLHEP_FIND_VERSION AND CLHEP_FIND_VERSION_EXACT)
            if ("CLHEP ${CLHEP_FIND_VERSION}" STREQUAL ${CLHEP_DISCOVERED_VERSION})
                set(CLHEP_VERSION_OK TRUE)
            else ()
                set(CLHEP_VERSION_OK FALSE)
            endif ()
        elseif (CLHEP_FIND_VERSION AND NOT CLHEP_FIND_VERSION_EXACT)
            if (CLHEP_MAX_VERSION)
                MACRO_ENSURE_VERSION_RANGE("${CLHEP_FIND_VERSION}" "${CLHEP_DISCOVERED_VERSION}" "${CLHEP_MAX_VERSION}" CLHEP_VERSION_OK)
            else ()
                MACRO_ENSURE_VERSION("${CLHEP_FIND_VERSION}" "${CLHEP_DISCOVERED_VERSION}" CLHEP_VERSION_OK)
            endif()
        else ()
            if (CLHEP_MAX_VERSION)
                MACRO_ENSURE_VERSION_RANGE("0.0.0.0" "${CLHEP_DISCOVERED_VERSION}" "${CLHEP_MAX_VERSION}" CLHEP_VERSION_OK)
            else()
                set(CLHEP_VERSION_OK TRUE)
            endif ()
        endif ()
        # END OF VERSIONING CHECKS

        # Use clhep-config to find libs and includes
        execute_process(COMMAND ${CLHEP_CONFIG_EXECUTABLE} --libs OUTPUT_VARIABLE CLHEP_CONFIG_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND ${CLHEP_CONFIG_EXECUTABLE} --include OUTPUT_VARIABLE CLHEP_CONFIG_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

        # Did we succeed?
        if (CLHEP_CONFIG_LIBS AND CLHEP_CONFIG_INCLUDE_DIR)
            #Resolve the includes and libraries
            RESOLVE_INCLUDES(CLHEP_INCLUDE_DIRS ${CLHEP_CONFIG_INCLUDE_DIR})
            RESOLVE_LIBRARIES(CLHEP_LIBRARIES ${CLHEP_CONFIG_LIBS})
        endif ()
    endif ()
endif ()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CLHEP "Failed to find CLHEP" CLHEP_VERSION_OK CLHEP_LIBRARIES CLHEP_INCLUDE_DIRS)
mark_as_advanced(CLHEP_CONFIG_EXECUTABLE)

