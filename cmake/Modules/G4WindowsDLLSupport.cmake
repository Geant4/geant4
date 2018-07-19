#.rst:
# G4WindowsDLLSupport
# -------------------
#
# Provides functions needed to support builds of DLLs both pre and post
# CMake 3.4. This version introduced the WINDOWS_EXPORT_ALL_SYMBOLS
# property, removing the need for the shipped `genwindef` program.
#
# With Geant4 10.4 and below having a minimum version requirement of
# CMake 3.3, we need to provide both so that patches can be applied
# to these versions without requiring clients to bump their CMake version.
#

#-----------------------------------------------------------------
# License and Disclaimer
#
# The  Geant4 software  is  copyright of the Copyright Holders  of
# the Geant4 Collaboration.  It is provided  under  the terms  and
# conditions of the Geant4 Software License,  included in the file
# LICENSE and available at  http://cern.ch/geant4/license .  These
# include a list of copyright holders.
#
# Neither the authors of this software system, nor their employing
# institutes,nor the agencies providing financial support for this
# work  make  any representation or  warranty, express or implied,
# regarding  this  software system or assume any liability for its
# use.  Please see the license in the file  LICENSE  and URL above
# for the full disclaimer and the limitation of liability.
#
# This  code  implementation is the result of  the  scientific and
# technical work of the GEANT4 collaboration.
# By using,  copying,  modifying or  distributing the software (or
# any work based  on the software)  you  agree  to acknowledge its
# use  in  resulting  scientific  publications,  and indicate your
# acceptance of all terms of the Geant4 Software license.
#
#-----------------------------------------------------------------

if(NOT __G4WINDOWSDLLSUPPORT_INCLUDED)
  set(__G4WINDOWSDLLSUPPORT_INCLUDED)
else()
  return()
endif()

# On WIN32, we need to build the genwindef application to create export
# def files for building DLLs.
# We only use it as a helper application at the moment so we exclude it from
# the ALL target.
if(WIN32)
  # Assume the sources are co-located
  get_filename_component(_genwindef_src_dir ${CMAKE_CURRENT_LIST_FILE} PATH)
  add_executable(genwindef EXCLUDE_FROM_ALL
    ${_genwindef_src_dir}/genwindef/genwindef.cpp
    ${_genwindef_src_dir}/genwindef/LibSymbolInfo.h
    ${_genwindef_src_dir}/genwindef/LibSymbolInfo.cpp)
endif()

#.rst:
# Macro for supporting DLL build on CMake < 3.4
macro(__geant4_add_dll_old)
  # It's an error to call this on non-WIN32 platforms
  if(NOT WIN32)
    message(FATAL_ERROR "__geant4_add_dll_old is not callable on non-Windows platforms!")
  endif()

  # - Takes same args as geant4_library_target for compatibility
  cmake_parse_arguments(__G4ADDOLDDLL
    ""
    "NAME" "SOURCES;GEANT4_LINK_LIBRARIES;LINK_LIBRARIES"
    ${ARGN}
    )

  # We have to generate the def export file from an archive library.
  # This is a temporary separate from a real archive library, and
  # even though it's static, we need to mark that it will have
  # DLL symbols via the G4LIB_BUILD_DLL macro
  add_library(_${__G4ADDOLDDLL_NAME}-archive STATIC EXCLUDE_FROM_ALL ${__G4ADDOLDDLL_SOURCES})
  set(_archive _${__G4ADDOLDDLL_NAME}-archive)
  target_compile_features(${_archive} PUBLIC ${GEANT4_TARGET_COMPILE_FEATURES})
  target_compile_definitions(${_archive} PUBLIC -DG4LIB_BUILD_DLL)

  # - Add the config specific compile definitions
  geant4_compile_definitions_config(${_archive})

  # - Create the .def file for this library
  # Use generator expressions to get path to per-mode lib and
  # older CMAKE_CFG_INTDIR variable to set name of per-mode def
  # file (Needed as generator expressions cannot be used in argument
  # to OUTPUT...
  add_custom_command(OUTPUT _${__G4ADDOLDDLL_NAME}-${CMAKE_CFG_INTDIR}.def
    COMMAND genwindef -o _${__G4ADDOLDDLL_NAME}-${CMAKE_CFG_INTDIR}.def -l ${__G4ADDOLDDLL_NAME} $<TARGET_FILE:${_archive}>
    DEPENDS ${_archive} genwindef)

  # - Now we can build the DLL
  # We create it from a dummy empty C++ file plus the def file.
  # Also set the public compile definition on it so that clients
  # will set correct macro automatically.
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/_${__G4ADDOLDDLL_NAME}.cpp
    "// empty _${__G4ADDOLDDLL_NAME}.cpp\n")

  add_library(${__G4ADDOLDDLL_NAME} SHARED _${__G4ADDOLDDLL_NAME}.cpp
    _${__G4ADDOLDDLL_NAME}-${CMAKE_CFG_INTDIR}.def)
  target_compile_definitions(${__G4ADDOLDDLL_NAME} PUBLIC -DG4LIB_BUILD_DLL)
  target_compile_features(${__G4ADDOLDDLL_NAME} PUBLIC ${GEANT4_TARGET_COMPILE_FEATURES})

  # - Link the DLL.
  # We link it to the archive, and the supplied libraries,
  # but then remove the archive from the LINK_INTERFACE.
  target_link_libraries(${__G4ADDOLDDLL_NAME}
    ${_archive}
    ${__G4ADDOLDDLL_GEANT4_LINK_LIBRARIES}
    ${__G4ADDOLDDLL_LINK_LIBRARIES})

  set_target_properties(${__G4ADDOLDDLL_NAME}
    PROPERTIES INTERFACE_LINK_LIBRARIES "${__G4ADDOLDDLL_GEANT4_LINK_LIBRARIES};${__G4ADDOLDDLL_LINK_LIBRARIES}")
endmacro()
