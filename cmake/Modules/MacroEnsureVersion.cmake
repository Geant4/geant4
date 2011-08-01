# This file defines the following macros for developers to use in ensuring
# that installed software is of the right version:
#
# MACRO_ENSURE_VERSION        - test that a version number is greater than
#                               or equal to some minimum
# MACRO_ENSURE_VERSION_RANGE - test that a version number is greater than
#                               or equal to some minimum and less than some
#                               maximum
# MACRO_ENSURE_VERSION2       - deprecated, do not use in new code
#
# MACRO_VERIFY_VERSION        - test that a version number satisfies a given
#                               set of requirements, such as an exact match,
#                               greater than or equal to a minimum, less than
#                               or equal to some maximum, or between a minimum
#                               and maximum.

# MACRO_ENSURE_VERSION
# This macro compares version numbers of the form "a.b.c.d" or "a.b.c" or "a.b"
# MACRO_ENSURE_VERSION( FOO_MIN_VERSION FOO_VERSION_FOUND FOO_VERSION_OK)
# will set FOO_VERSION_OK to true if FOO_VERSION_FOUND >= FOO_MIN_VERSION
# Leading and trailing text is ok, e.g.
# MACRO_ENSURE_VERSION( "2.5.31" "flex 2.5.4a" VERSION_OK)
# which means 2.5.31 is required and "flex 2.5.4a" is what was found on the system

# Copyright (c) 2006, David Faure, <faure@kde.org>
# Copyright (c) 2007, Will Stephenson <wstephenson@kde.org>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.


# MACRO_VERIFY_VERSION
# This macro compares version numbers of the form "a.b.c.d" or "a.b.c" or "a.b"
#
# MACRO_VERIFY_VERSION(VERSION_OK 
#   MIN_VERSION "a.b.c.d"
#   MAX_VERSION "p.q.r"
#   FOUND_VERSION "x.y.z"
#   [FIND_EXACT True|False])
#
# The macro uses the following logic to verify version checking
#
# 1) MIN_VERSION a.b.c, FOUND_VERSION x.y.z supplied, FIND_EXACT set:
#    Check and require that x.y.z == a.b.c
#    MAX_VERSION ignored.
#
# 2) MIN_VERSION a.b.c, FOUND_VERSION x.y.z supplied
#   i)  If MAX_VERSION not set, check and require that x.y.z >= a.b.c
#   ii) If MAX_VERSION p.q.r set, check and require a.b.c <= x.y.z < p.q.r
#
# 3) MIN_VERSION a.b.c NOT supplied, FOUND_VERSION x.y.z supplied
#   i)  If MAX_VERSION not set, version assumed to be verified.
#   ii) If MAX_VERSION p.q.r set, check and require x.y.z < p.q.r
#
# Copyright (c) 2009, Ben Morgan, <Ben.Morgan@warwick.ac.uk>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.




# MACRO_ENSURE_VERSION_RANGE
# This macro ensures that a version number of the form
# "a.b.c.d" or "a.b.c" or "a.b" falls within a range defined by
# min_version <= found_version < max_version.
# If this expression holds, FOO_VERSION_OK will be set TRUE
#
# Example: MACRO_ENSURE_VERSION_RANGE3( "0.1.0" ${FOOCODE_VERSION} "0.7.0" FOO_VERSION_OK )
#
# This macro will break silently if any of a,b,c or d are greater than 100.
#
# Copyright (c) 2007, Will Stephenson <wstephenson@kde.org>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

# NORMALIZE_VERSION
# Helper macro to convert version numbers of the form "a.b.c.d"
# to an integer equal to 10^6 * a + 10^4 * b + 10^2 * c + d
#
# This macro will break silently if any of a,b,c,d are greater than 100.
#
# Copyright (c) 2006, David Faure, <faure@kde.org>
# Copyright (c) 2007, Will Stephenson <wstephenson@kde.org>
#
# Modifications Copyright (c) 2008, Ben Morgan <Ben.Morgan@warwick.ac.uk>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

# CHECK_RANGE_INCLUSIVE_LOWER
# Helper macro to check whether x <= y < z
#
# Copyright (c) 2007, Will Stephenson <wstephenson@kde.org>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

# MACRO_

MACRO(NORMALIZE_VERSION _requested_version _normalized_version)
    # Do we have a.b.c.d in _requested_version?
    STRING(REGEX MATCH "[^0-9]*[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9].*" _fourPartMatch "${_requested_version}")
    STRING(REGEX MATCH "[^0-9]*[0-9]+\\.[0-9]+\\.[0-9]+.*" _threePartMatch "${_requested_version}")

    if (_fourPartMatch)
        # We certainly have a.b.c.d so split out the components as needed
        STRING(REGEX REPLACE "[^0-9]*([0-9]+)\\.[0-9]+\\.[0-9]+.*" "\\1" _major_vers "${_requested_version}")
        STRING(REGEX REPLACE "[^0-9]*[0-9]+\\.([0-9]+)\\.[0-9]+.*" "\\1" _minor_vers "${_requested_version}")
        STRING(REGEX REPLACE "[^0-9]*[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" _patch_vers "${_requested_version}")
        STRING(REGEX REPLACE "[^0-9]*[0-9]+\\.[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" _tweak_vers "${_requested_version}")

    elseif (_threePartMatch)
        # parse the parts of the version string
        STRING(REGEX REPLACE "[^0-9]*([0-9]+)\\.[0-9]+\\.[0-9]+.*" "\\1" _major_vers "${_requested_version}")
        STRING(REGEX REPLACE "[^0-9]*[0-9]+\\.([0-9]+)\\.[0-9]+.*" "\\1" _minor_vers "${_requested_version}")
        STRING(REGEX REPLACE "[^0-9]*[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" _patch_vers "${_requested_version}")
        set(_tweak_vers "0")
    else (_threePartMatch)
        STRING(REGEX REPLACE "[^0-9]*([0-9]+)\\.[0-9]+" "\\1" _major_vers "${_requested_version}")
        STRING(REGEX REPLACE "[^0-9]*[0-9]+\\.([0-9]+)" "\\1" _minor_vers "${_requested_version}")
        set(_patch_vers "0")
        set(_tweak_vers "0")
    endif (_fourPartMatch)

    # compute an overall version number which can be compared at once
    MATH(EXPR ${_normalized_version} "${_major_vers}*1000000 + ${_minor_vers}*10000 + ${_patch_vers}*100 + ${_tweak_vers}")
ENDMACRO(NORMALIZE_VERSION)

MACRO(MACRO_CHECK_RANGE_INCLUSIVE_LOWER _lower_limit _value _upper_limit _ok)
   if (${_value} LESS ${_lower_limit})
      set( ${_ok} FALSE )
  elseif (${_value} EQUAL ${_lower_limit})
      set( ${_ok} TRUE )
  elseif (${_value} EQUAL ${_upper_limit})
      set( ${_ok} FALSE )
  elseif (${_value} GREATER ${_upper_limit})
      set( ${_ok} FALSE )
  else (${_value} LESS ${_lower_limit})
      set( ${_ok} TRUE )
  endif (${_value} LESS ${_lower_limit})
ENDMACRO(MACRO_CHECK_RANGE_INCLUSIVE_LOWER)

MACRO(MACRO_ENSURE_VERSION requested_version found_version var_too_old)
    NORMALIZE_VERSION( ${requested_version} req_vers_num )
    NORMALIZE_VERSION( ${found_version} found_vers_num )

    if (found_vers_num LESS req_vers_num)
        set( ${var_too_old} FALSE )
    else (found_vers_num LESS req_vers_num)
        set( ${var_too_old} TRUE )
    endif (found_vers_num LESS req_vers_num)

ENDMACRO(MACRO_ENSURE_VERSION)

MACRO(MACRO_ENSURE_VERSION2 requested_version2 found_version2 var_too_old2)
    MACRO_ENSURE_VERSION( ${requested_version2} ${found_version2} ${var_too_old2})
ENDMACRO(MACRO_ENSURE_VERSION2)

MACRO(MACRO_ENSURE_VERSION_RANGE min_version found_version max_version var_ok)
    NORMALIZE_VERSION( ${min_version} req_vers_num )
    NORMALIZE_VERSION( ${found_version} found_vers_num )
    NORMALIZE_VERSION( ${max_version} max_vers_num )

    MACRO_CHECK_RANGE_INCLUSIVE_LOWER( ${req_vers_num} ${found_vers_num} ${max_vers_num} ${var_ok})
ENDMACRO(MACRO_ENSURE_VERSION_RANGE)


include(CMakeMacroParseArguments)
MACRO(MACRO_VERIFY_VERSION)
    CMAKE_PARSE_ARGUMENTS(MVV
        ""
        "MIN_VERSION;MAX_VERSION;FOUND_VERSION;FIND_EXACT" 
        ""
        ${ARGN})

    set(MVV_VAR_OK ${ARGV0})
    message("var to set: ${ARGV0}")
    message("min_version = ${MVV_MIN_VERSION}")
    message("max_version = ${MVV_MAX_VERSION}")
    message("found_version = ${MVV_FOUND_VERSION}")
    message("find_exact = ${MVV_FIND_EXACT}")

    # If MIN_VERSION NOT supplied, must set it to absolute minimum
    if (NOT DEFINED ${MVV_MIN_VERSION})
        set(MVV_MIN_VERSION "0.0.0.0")
    endif (NOT DEFINED ${MVV_MIN_VERSION})

    if (MVV_FOUND_VERSION AND MVV_FIND_EXACT)
        if (${MVV_MIN_VERSION} MATCHES ${MVV_FOUND_VERSION})
            set(${MVV_VAR_OK} TRUE)
        else (${MVV_MIN_VERSION} MATCHES ${MVV_FOUND_VERSION})
            set(${MVV_VAR_OK} FALSE)
        endif (${MVV_MIN_VERSION} MATCHES ${MVV_FOUND_VERSION})
    elseif (MVV_FOUND_VERSION AND NOT MVV_FIND_EXACT)
        if (MVV_MAX_VERSION)
            MACRO_ENSURE_VERSION_RANGE(${MVV_MIN_VERSION}
                ${MVV_FOUND_VERSION}
                ${MVV_MAX_VERSION}
                ${MVV_VAR_OK})
        else (MVV_MAX_VERSION)
            MACRO_ENSURE_VERSION(${MVV_MIN_VERSION} 
                ${MVV_FOUND_VERSION} 
                ${MVV_VAR_OK})
        endif (MVV_MAX_VERSION)
    else (MVV_FOUND_VERSION AND MVV_FIND_EXACT)
        if (MVV_MAX_VERSION)
            MACRO_ENSURE_VERSION_RANGE("0.0.0.0"
                ${MVV_FOUND_VERSION}
                ${MVV_MAX_VERSION}
                ${MVV_VAR_OK})
        else (MVV_MAX_VERSION)
            set(${MVV_VAR_OK} TRUE)
        endif (MVV_MAX_VERSION)
    endif (MVV_FOUND_VERSION AND MVV_FIND_EXACT)
ENDMACRO(MACRO_VERIFY_VERSION)

