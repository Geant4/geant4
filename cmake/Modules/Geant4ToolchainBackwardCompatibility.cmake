# - Script for setting up backward compatibility with GNU make user toolchain
#
# The GNU make based buildsystem for Geant4 provides a toolchain for user
# building simple Geant4 applications. The build and install of Geant4 
# provides a customized set of non-standard install paths with use of the
# toolchain dependent on environment variables pointing to the install
# paths.
#
# This script processes information on the CMake install paths, system and
# compiler to determine the following variables for backward compatibility
#
#     GEANT4_SYSTEM  : Old style system name, e.g. 'Linux', 'Darwin' or 'WIN32'
#
#     GEANT4_COMPILER: Old system compiler id, e.g. 'g++', 'VC'. 
#
#     G4INSTALL      : Location of config subdirectory containing GNU make
#                      fragment files for toolchain.
#
#     G4INCLUDE      : Old style path to location of Geant4 headers
#
#     G4LIB          : Old style library directory path. Rather than
#                      containing the actual libraries, it is expected to
#                      contain subdirectories named
#                      GEANT4_SYSTEM-GEANT4_COMPILER
#
#
# These variables are used in CMake configuration files which generate shell
# scripts (C and Bourne flavour) the user can source to set up their
# environment. These replace the old style 'env.(c)sh' scripts to allow
# users to work with the new CMake built libraries transparently if
# their application relies on the old style toolchain.
#
# $Id: Geant4ToolchainBackwardCompatibility.cmake,v 1.4 2010-12-01 17:13:20 bmorgan Exp $
# GEANT4 Tag $Name: not supported by cvs2svn $
#


#------------------------------------------------------------------------------
# Determine the backward compatible system name
#
if(NOT WIN32)
    set(GEANT4_SYSTEM ${CMAKE_SYSTEM_NAME})
else()
    set(GEANT4_SYSTEM "WIN32")
endif()

message(STATUS "Geant4 backwards compatible system name: ${GEANT4_SYSTEM}")

#------------------------------------------------------------------------------
# Determine the backward compatible compiler name
#
if(CMAKE_COMPILER_IS_GNUCXX)
    set(GEANT4_COMPILER "g++")
elseif(MSVC)
    set(GEANT4_COMPILER "VC")
elseif(CMAKE_CXX_COMPILER MATCHES "icpc.*")
    set(GEANT4_COMPILER "icc")
else()
    set(GEANT4_COMPILER "UNSUPPORTED")
endif()

message(STATUS "Geant4 backwards compatible compiler name: ${GEANT4_COMPILER}")


#------------------------------------------------------------------------------
# Construct backward compatible variables
#
set(G4SYSTEM  "${GEANT4_SYSTEM}-${GEANT4_COMPILER}")
set(G4INSTALL ${GEANT4_DATADIR}/geant4-${geant4_VERSION})
set(G4INCLUDE ${GEANT4_INCLUDEDIR}/geant4)
set(G4LIB     ${GEANT4_LIBDIR}/geant4-${geant4_VERSION})

message(STATUS "Geant4 backwards compatible variable G4SYSTEM : ${G4SYSTEM}")
message(STATUS "Geant4 backwards compatible variable G4INSTALL: ${G4INSTALL}")
message(STATUS "Geant4 backwards compatible variable G4INCLUDE: ${G4INCLUDE}")
message(STATUS "Geant4 backwards compatible variable G4LIB    : ${G4LIB}")


#------------------------------------------------------------------------------
# CLHEP setup
# Unless CLHEP is setup very oddly, the base directory should be sufficient
execute_process(COMMAND ${CLHEP_CONFIG_EXECUTABLE} --prefix 
    OUTPUT_VARIABLE CLHEP_BASE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)



#------------------------------------------------------------------------------
# Optional variables which are set only if certain features are enabled
#
if(GEANT4_USE_QT)
    set(G4_HAVE_QT "yes")
else()
    set(G4_HAVE_QT "no")
endif()

if(UNIX)
    set(G4_HAVE_TCSH "yes")
else()
    set(G4_HAVE_TCSH "no")
endif()



