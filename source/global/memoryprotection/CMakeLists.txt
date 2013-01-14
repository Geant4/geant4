#-----------------------------------------------------------------------
# memoryprotection library for Geant4
#
# memoryprotection is built as a seperate project to allow easy movement 
# and/or integration as development progresses.
# The location where the libraries and headers are to be installed can
# be adjusted with the following variables (assumed relative to
# CMAKE_INSTALL_PREFIX):
#
# MEMORYPROTECTION_INSTALL_BINDIR       The installation directory for 
#                                       executables
# MEMORYPROTECTION_INSTALL_INCLUDEDIR   The installation directory for 
#                                       headers
# MEMORYPROTECTION_INSTALL_LIBDIR       The installation directory for 
#                                       libraries
# 
# MEMORYPROTECTION_INSTALL_EXPORT_NAME  The EXPORT option value for 
#                                       install(TARGETS) calls.
#

#-----------------------------------------------------------------------
# Set up the project
#
cmake_minimum_required(VERSION 2.6.4 FATAL_ERROR)
project(memoryprotection C)

#-----------------------------------------------------------------------
# Fail if we don't have GNU C compiler
#
if(NOT CMAKE_C_COMPILER_ID STREQUAL GNU)
  message(FATAL_ERROR "memoryprotection requires the GNU C compiler")
endif()

#-----------------------------------------------------------------------
# Need pthreads
#
find_package(Threads)

#-----------------------------------------------------------------------
# Force the compiler flags to always use -g -O0
# Handle both build and configuration types
set(CMAKE_C_FLAGS "-g -O0")

string(TOUPPER "${CMAKE_BUILD_TYPE}" MEMORYPROTECTION_BUILD_TYPE)
set(CMAKE_C_FLAGS_${MEMORYPROTECTION_BUILD_TYPE} "")

foreach(_type ${CMAKE_CONFIGURATION_TYPES})
  string(TOUPPER "${_type}" _upper_type)
  set(CMAKE_C_FLAGS_${MEMORYPROTECTION_BUILD_TYPE} "")
endforeach()

#-----------------------------------------------------------------------
# If MEMORYPROTECTION_INSTALL_<COMPONENT>DIR are not set, provide defaults
#

if(NOT MEMORYPROTECTION_INSTALL_BINDIR)
  set(MEMORYPROTECTION_INSTALL_BINDIR "bin")
endif()

if(NOT MEMORYPROTECTION_INSTALL_INCLUDEDIR)
  set(MEMORYPROTECTION_INSTALL_INCLUDEDIR "include")
endif()

if(NOT MEMORYPROTECTION_INSTALL_LIBDIR)
  set(MEMORYPROTECTION_INSTALL_LIBDIR "lib")
endif()

if(NOT MEMORYPROTECTION_INSTALL_EXPORT_NAME)
  set(MEMORYPROTECTION_INSTALL_EXPORT_NAME "memoryprotectionlibrarydepends")
endif()


#-----------------------------------------------------------------------
# Build the MemoryProtection library, always shared
# NB we include the header in the source list so that IDEs will include
# it in the solutions
add_library(ProtectedMemory SHARED ProtectedMemory.h ProtectedMemory.c)

#-----------------------------------------------------------------------
# Build the tracer program
#
add_executable(tracer tracer.c filter.h filter.c MemorySegment.h UDSignals.h)
target_link_libraries(tracer ${CMAKE_THREAD_LIBS_INIT})

#-----------------------------------------------------------------------
# Install commands
#
install(TARGETS ProtectedMemory tracer
  EXPORT ${MEMORYPROTECTION_INSTALL_EXPORT_NAME}
  RUNTIME DESTINATION ${MEMORYPROTECTION_INSTALL_BINDIR}
  LIBRARY DESTINATION ${MEMORYPROTECTION_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${MEMORYPROTECTION_INSTALL_LIBDIR})

install(FILES ProtectedMemory.h
  DESTINATION "${MEMORYPROTECTION_INSTALL_INCLUDEDIR}/memoryprotection")

