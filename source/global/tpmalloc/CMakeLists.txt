#-----------------------------------------------------------------------
# tpmalloc library for Geant4
#
# tpmalloc is built as a seperate project to allow easy movement and/or
# integration as development progresses.
# The location where the libraries and headers are to be installed can
# be adjusted with the following variables (assumed relative to
# CMAKE_INSTALL_PREFIX):
#
# TPMALLOC_INSTALL_BINDIR        The installation directory for executables
# TPMALLOC_INSTALL_INCLUDEDIR    The installation directory for headers
# TPMALLOC_INSTALL_LIBDIR        The installation directory for libraries
# 
# TPMALLOC_INSTALL_EXPORT_NAME   The EXPORT option value for 
#                                install(TARGETS) calls.
#

#-----------------------------------------------------------------------
# Set up the project
#
cmake_minimum_required(VERSION 2.6.4 FATAL_ERROR)
project(tpmalloc C)

#-----------------------------------------------------------------------
# Fail if we don't have GNU C compiler
#
if(NOT CMAKE_C_COMPILER_ID STREQUAL GNU)
  message(FATAL_ERROR "tpmalloc requires the GNU C compiler")
endif()

#-----------------------------------------------------------------------
# Force the compiler flags to always use -O3
# Handle both build and configuration types
set(CMAKE_C_FLAGS "-O3")

string(TOUPPER "${CMAKE_BUILD_TYPE}" TPMALLOC_BUILD_TYPE)
set(CMAKE_C_FLAGS_${TPMALLOC_BUILD_TYPE} "")

foreach(_type ${CMAKE_CONFIGURATION_TYPES})
  string(TOUPPER "${_type}" _upper_type)
  set(CMAKE_C_FLAGS_${TPMALLOC_BUILD_TYPE} "")
endforeach()

#-----------------------------------------------------------------------
# If TPMALLOC_INSTALL_<COMPONENT>DIR are not set, provide defaults
#

if(NOT TPMALLOC_INSTALL_BINDIR)
  set(TPMALLOC_INSTALL_BINDIR "bin")
endif()

if(NOT TPMALLOC_INSTALL_INCLUDEDIR)
  set(TPMALLOC_INSTALL_INCLUDEDIR "include")
endif()

if(NOT TPMALLOC_INSTALL_LIBDIR)
  set(TPMALLOC_INSTALL_LIBDIR "lib")
endif()

if(NOT TPMALLOC_INSTALL_EXPORT_NAME)
  set(TPMALLOC_INSTALL_EXPORT_NAME "tpmalloclibrarydepends")
endif()


#-----------------------------------------------------------------------
# Build the tpmallocstub and mymalloc libs, always shared
#
add_library(tpmallocstub SHARED tpmallocstub.h tpmallocstub.c)
add_library(mymalloc SHARED hjmalloc.c mymalloc.h mymalloc.c)

#-----------------------------------------------------------------------
# Install commands
#
install(TARGETS tpmallocstub mymalloc
  EXPORT ${TPMALLOC_INSTALL_EXPORT_NAME}
  RUNTIME DESTINATION ${TPMALLOC_INSTALL_BINDIR}
  LIBRARY DESTINATION ${TPMALLOC_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${TPMALLOC_INSTALL_LIBDIR})

install(FILES tpmallocstub.h mymalloc.h
  DESTINATION "${TPMALLOC_INSTALL_INCLUDEDIR}/tpmalloc")

