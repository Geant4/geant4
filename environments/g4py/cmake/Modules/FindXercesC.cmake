# - Find Xerces-C
# This module tries to find the Xerces-C library and headers.
# Once done this will define
#
#   XERCESC_FOUND - system has Xerces-C headers and libraries
#   XERCESC_INCLUDE_DIRS - the include directories needed for Xerces-C
#   XERCESC_LIBRARIES - the libraries needed to use Xerces-C
#
# Variables used by this module, which can change the default behaviour and
# need to be set before calling find_package:
#
#   XERCESC_ROOT_DIR            Root directory to Xerces-C installation. Will
#                               be used ahead of CMake default path.
#
# The following advanced variables may be used if the module has difficulty
# locating Xerces-C or you need fine control over what is used.
#
#   XERCESC_INCLUDE_DIR
#
#   XERCESC_LIBRARY
#
# Copyright (c) 2009, Ben Morgan, <Ben.Morgan@warwick.ac.uk>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.



# Look for the header - preferentially searching below XERCESC_ROOT_DIR
find_path(
    XERCESC_INCLUDE_DIR 
    NAMES xercesc/util/XercesVersion.hpp
    PATHS ${XERCESC_ROOT_DIR}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
)

# If we didn't find it there, fall back to some standard search paths
find_path(
    XERCESC_INCLUDE_DIR
    NAMES xercesc/util/XercesVersion.hpp
)

# Look for the library, preferentially searching below XERCESC_ROOT_DIR
find_library(
    XERCESC_LIBRARY
    NAMES xerces-c xerces-c_3
    PATHS ${XERCESC_ROOT_DIR}
    PATH_SUFFIXES lib64 lib32 lib
    NO_DEFAULT_PATH
)

find_library(
    XERCESC_LIBRARY
    NAMES xerces-c xerces-c_3
)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    XercesC
    DEFAULT_MSG
    XERCESC_LIBRARY
    XERCESC_INCLUDE_DIR
)

if (XERCESC_FOUND)
    set(XERCESC_LIBRARIES ${XERCESC_LIBRARY})
    set(XERCESC_INCLUDE_DIRS ${XERCESC_INCLUDE_DIR})
else (XERCESC_FOUND)
    set(XERCESC_LIBRARIES)
    set(XERCESC_INCLUDE_DIRS)
endif (XERCESC_FOUND)


mark_as_advanced(
    XERCESC_LIBRARY 
    XERCESC_INCLUDE_DIR
)
