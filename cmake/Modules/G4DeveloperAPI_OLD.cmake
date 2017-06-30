#.rst:
# G4DeveloperAPI_OLD
# ------------------
#
# Old style API for declaring Geant4 source code modules and
# compiling, linking and installing binary libraries composed
# from one or more source code modules.
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


# - Include guard
if(NOT __G4DEVELOPERAPI_OLD_INCLUDED)
  set(__G4DEVELOPERAPI_OLD_INCLUDED 1)
else()
  return()
endif()

# - Needed CMake modules
include(CMakeParseArguments)

#-----------------------------------------------------------------------
#.rst:
# Source Code Module Declaration Functions
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Geant4's source code is organised into categories composed of one or
# more "modules" of source code. Each module has the following layout
# of headers and sources:
#
# - include/
#   - Public headers
# - src/
#   - Implementations
#
# The functions here declare a module ready for it to be integrated
# by the library build functions described below.

#-----------------------------------------------------------------------
# macro geant4_define_module(NAME <name>
#                            HEADERS <header1> <header2> ... <headerN>
#                            SOURCES <source1> <source2> ... <sourceN>
#                            GRANULAR_DEPENDENCIES <dep1> ... <depN>
#                            GLOBAL_DEPENDENCIES <dep1> ... <depN>
#                            LINK_LIBRARIES <lib1> ... <lib2>)
#       Define a Geant4 Module's sources and what internal and external
#       libraries it links to.
#
# GEANT4_DEFINE_MODULE will take the name of the module, a list of header
# files, a list of source files and dependencies:
#
# GEANT4_DEFINE_MODULE(NAME name
#                      HEADERS hdr1 hdr2 ...
#                      SOURCES src1 src2 ...
#                      GRANULAR_DEPENDENCIES dep1 dep2 ...
#                      GLOBAL_DEPENDENCIES dep1 dep2 ...
#                      LINK_LIBRARIES lib1 lib2 ...)
#
# It assumes that it will be called from a CMake located in the module
# root directory (i.e. the directory containing the include and src
# subdirectories). It uses this location to define absolute paths to the
# headers and sources, and will ignore any absolute paths passed in HEADERS
# or SOURCES so that generated files can be passed in.
#
# The macro defines the variables
#
# G4MODULENAME = name
# ${G4MODULENAME}_HEADERS = List of absolute paths to files given in HEADERS
# ${G4MODULENAME}_SOURCES = List of absolute paths to files given in SOURCES
# ${G4MODULENAME}_GRANULAR_DEPS  = List of granular libraries on which this
#                                  module depends
# ${G4MODULENAME}_GLOBAL_DEPS    = List of global libraries on which this module
#                                  depends
# ${G4MODULENAME}_LINK_LIBRARIES = List of external libraries to be linked to
#                                  this module
#
# It will also add the module include directory to the list of directories
# using include_directories
#
macro(geant4_define_module)
  cmake_parse_arguments(G4DEFMOD
    ""
    "NAME"
    "HEADERS;SOURCES;GRANULAR_DEPENDENCIES;GLOBAL_DEPENDENCIES;LINK_LIBRARIES"
    ${ARGN}
    )

  set(G4MODULENAME ${G4DEFMOD_NAME})

  get_filename_component(${G4MODULENAME}_BASEDIR ${CMAKE_CURRENT_LIST_FILE} PATH)
  set(${G4MODULENAME}_SRCDIR ${${G4MODULENAME}_BASEDIR}/src)
  set(${G4MODULENAME}_INCDIR ${${G4MODULENAME}_BASEDIR}/include)

  # We now create absolute paths to the headers and sources,
  # ignoring any already absolute paths
  foreach(_HDR ${G4DEFMOD_HEADERS})
    if(IS_ABSOLUTE ${_HDR})
      list(APPEND ${G4MODULENAME}_HEADERS ${_HDR})
    else()
      list(APPEND ${G4MODULENAME}_HEADERS ${${G4MODULENAME}_INCDIR}/${_HDR})
    endif()
  endforeach()

  foreach(_SRC ${G4DEFMOD_SOURCES})
    if(IS_ABSOLUTE ${_SRC})
      list(APPEND ${G4MODULENAME}_SOURCES ${_SRC})
    else()
      list(APPEND ${G4MODULENAME}_SOURCES ${${G4MODULENAME}_SRCDIR}/${_SRC})
    endif()
  endforeach()

  foreach(_LIB ${G4DEFMOD_GRANULAR_DEPENDENCIES})
    list(APPEND ${G4MODULENAME}_GRANULAR_DEPENDENCIES ${_LIB})
  endforeach()

  foreach(_LIB ${G4DEFMOD_GLOBAL_DEPENDENCIES})
    list(APPEND ${G4MODULENAME}_GLOBAL_DEPENDENCIES ${_LIB})
  endforeach()

  foreach(_LIB ${G4DEFMOD_LINK_LIBRARIES})
    list(APPEND ${G4MODULENAME}_LINK_LIBRARIES ${_LIB})
  endforeach()

  include_directories(${${G4MODULENAME}_INCDIR})
endmacro()


#-----------------------------------------------------------------------
# macro geant4_add_compile_definitions(SOURCES <source1> ... <sourceN>
#                                      COMPILE_DEFINITIONS <def1> ... <defN>
#                                      )
#       Add extra compile definitions to a specific list of sources
#       in the current module. Macroized to handle the need to specify absolute paths.
#       and *must* be called at the same level as geant4_define_module
#
macro(geant4_add_compile_definitions)
  cmake_parse_arguments(G4ADDDEF
    ""
    ""
    "SOURCES;COMPILE_DEFINITIONS"
    ${ARGN}
    )

  # We assume that the sources have been added at the level of a
  # a sources.cmake, so are inside the src subdir of the sources.cmake
  get_filename_component(_ACD_BASE_PATH ${CMAKE_CURRENT_LIST_FILE} PATH)

  # Now for each file, add the definitions
  foreach(_acd_source ${G4ADDDEF_SOURCES})
    # Extract any existing compile definitions
    get_source_file_property(
      _acd_existing_properties
      ${_ACD_BASE_PATH}/src/${_acd_source}
      COMPILE_DEFINITIONS)

    if(_acd_existing_properties)
      set(_acd_new_defs ${_acd_existing_properties}
        ${G4ADDDEF_COMPILE_DEFINITIONS})
    else()
      set(_acd_new_defs ${G4ADDDEF_COMPILE_DEFINITIONS})
    endif()

    # quote compile defs because this must epand to space separated list
    set_source_files_properties(${_ACD_BASE_PATH}/src/${_acd_source}
      PROPERTIES COMPILE_DEFINITIONS "${_acd_new_defs}")
  endforeach()
endmacro()


#-----------------------------------------------------------------------
#.rst:
# Library Declaration, Build and Install API
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# - Define useful macros for building and installing Geant4 library targets
#
# This file defines the following macros for Geant4 developers needing to
# add shared and static library targets.
#
# GEANT4_LIBRARY_TARGET        - define standard Geant4 library targets
#
# The macro will take the name of the library and its sources, defining
# static and shared targets depending on the value of BUILD_SHARED_LIBS and
# BUILD_STATIC_LIBS. Install targets are also created.
#
# A custom compile definition "GEANT4_DEVELOPER_<CONFIG>" is set on
# each target using the target property COMPILE_DEFINITIONS_<CONFIG>
# target property


#-----------------------------------------------------------------------
# function geant4_compile_definitions_config(<target>)
#          Set a custom compile definition for a Geant4 target on a
#          per configuration basis:
#            For mode <CONFIG>, define GEANT4_DEVELOPER_<CONFIG>
#
function(geant4_compile_definitions_config _target)
  if(NOT TARGET ${_target})
    message(FATAL_ERROR "geant4_compile_definitions_config passed target '${_target}' which is not a valid CMake target")
  endif()

  # CMake 3 prefers $<CONFIG> in generator expressions
  set(__default_config_generator "\$<CONFIG>")

  set_property(TARGET ${_target}
    APPEND PROPERTY COMPILE_DEFINITIONS
      GEANT4_DEVELOPER_${__default_config_generator}
    )
endfunction()


#-----------------------------------------------------------------------
# - GEANT4_LIBRARY_TARGET
# General build and install of a Geant4 library target
#
macro(geant4_library_target)
  cmake_parse_arguments(G4LIBTARGET
    ""
    "NAME" "SOURCES;GEANT4_LINK_LIBRARIES;LINK_LIBRARIES"
    ${ARGN}
    )

  if(BUILD_SHARED_LIBS)
    # Add the shared library target and link its dependencies
    # WIN32 first
    if(WIN32)
      # We have to generate the def export file from an archive library.
      # This is a temporary separate from a real archive library, and
      # even though it's static, we need to mark that it will have
      # DLL symbols via the G4LIB_BUILD_DLL macro
      add_library(_${G4LIBTARGET_NAME}-archive STATIC EXCLUDE_FROM_ALL ${G4LIBTARGET_SOURCES})
      set(_archive _${G4LIBTARGET_NAME}-archive)
      target_compile_features(${_archive} PUBLIC ${GEANT4_TARGET_COMPILE_FEATURES})
      target_compile_definitions(${_archive} PUBLIC -DG4LIB_BUILD_DLL)

      # - Add the config specific compile definitions
      geant4_compile_definitions_config(${_archive})

      # - Create the .def file for this library
      # Use generator expressions to get path to per-mode lib and
      # older CMAKE_CFG_INTDIR variable to set name of per-mode def
      # file (Needed as generator expressions cannot be used in argument
      # to OUTPUT...
      add_custom_command(OUTPUT _${G4LIBTARGET_NAME}-${CMAKE_CFG_INTDIR}.def
        COMMAND genwindef -o _${G4LIBTARGET_NAME}-${CMAKE_CFG_INTDIR}.def -l ${G4LIBTARGET_NAME} $<TARGET_FILE:${_archive}>
        DEPENDS ${_archive} genwindef)

      # - Now we can build the DLL
      # We create it from a dummy empty C++ file plus the def file.
      # Also set the public compile definition on it so that clients
      # will set correct macro automatically.
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/_${G4LIBTARGET_NAME}.cpp
        "// empty _${G4LIBTARGET_NAME}.cpp\n")

      add_library(${G4LIBTARGET_NAME} SHARED _${G4LIBTARGET_NAME}.cpp
        _${G4LIBTARGET_NAME}-${CMAKE_CFG_INTDIR}.def)
      target_compile_definitions(${G4LIBTARGET_NAME} PUBLIC -DG4LIB_BUILD_DLL)
      target_compile_features(${G4LIBTARGET_NAME} PUBLIC ${GEANT4_TARGET_COMPILE_FEATURES})

      # - Link the DLL.
      # We link it to the archive, and the supplied libraries,
      # but then remove the archive from the LINK_INTERFACE.
      target_link_libraries(${G4LIBTARGET_NAME}
        ${_archive}
        ${G4LIBTARGET_GEANT4_LINK_LIBRARIES}
        ${G4LIBTARGET_LINK_LIBRARIES})

      set_target_properties(${G4LIBTARGET_NAME}
        PROPERTIES INTERFACE_LINK_LIBRARIES "${G4LIBTARGET_GEANT4_LINK_LIBRARIES};${G4LIBTARGET_LINK_LIBRARIES}")

    else()
      # - We build a Shared library in the usual fashion...
      add_library(${G4LIBTARGET_NAME} SHARED ${G4LIBTARGET_SOURCES})
      geant4_compile_definitions_config(${G4LIBTARGET_NAME})
      target_compile_features(${G4LIBTARGET_NAME} PUBLIC ${GEANT4_TARGET_COMPILE_FEATURES})
      target_link_libraries(${G4LIBTARGET_NAME}
        ${G4LIBTARGET_GEANT4_LINK_LIBRARIES}
        ${G4LIBTARGET_LINK_LIBRARIES})
    endif()

    # This property is set to prevent concurrent builds of static and
    # shared libs removing each others files.
    set_target_properties(${G4LIBTARGET_NAME}
      PROPERTIES CLEAN_DIRECT_OUTPUT 1)

    # Always use '@rpath' in install names of libraries. This is the
    # most flexible mechanism for
    set_target_properties(${G4LIBTARGET_NAME}
      PROPERTIES MACOSX_RPATH 1
      )

    # Install the library - note the use of RUNTIME, LIBRARY and ARCHIVE
    # this helps with later DLL builds.
    # Export to standard depends file for later install
    # NEEDS WORK TO REMOVE HARDCODED LIB/BIN DIR
    install(TARGETS ${G4LIBTARGET_NAME}
      EXPORT Geant4LibraryDepends
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

    # Append the library target to a global property so that build tree
    # export of library dependencies can pick up all targets
    set_property(GLOBAL APPEND
      PROPERTY GEANT4_EXPORTED_TARGETS ${G4LIBTARGET_NAME})
  endif()

  #
  # As above, but for static rather than shared library
  if(BUILD_STATIC_LIBS)
    # We have to distinguish the static from shared lib, so use -static in
    # name. Link its dependencies, and ensure we actually link to the
    # -static targets (We should strictly do this for the external
    # libraries as well if we want a pure static build).
    add_library(${G4LIBTARGET_NAME}-static STATIC ${G4LIBTARGET_SOURCES})
    geant4_compile_definitions_config(${G4LIBTARGET_NAME}-static)
    target_compile_features(${G4LIBTARGET_NAME}-static PUBLIC ${GEANT4_TARGET_COMPILE_FEATURES})

    set(G4LIBTARGET_GEANT4_LINK_LIBRARIES_STATIC )
    foreach(_tgt ${G4LIBTARGET_GEANT4_LINK_LIBRARIES})
      list(APPEND G4LIBTARGET_GEANT4_LINK_LIBRARIES_STATIC ${_tgt}-static)
    endforeach()

    # If we are building both types of library and builtin clhep etc,
    # we want to link shared->shared and static->static.
    # Because externals like clhep appear in G4LIBTARGET_LINK_LIBRARIES,
    # filter this list to replace shared builtins with their static variant
    string(REGEX REPLACE
      "(G4clhep|G4expat|G4zlib|G4geomUSolids)(;|$)" "\\1-static\\2"
      G4LIBTARGET_LINK_LIBRARIES_STATIC
      "${G4LIBTARGET_LINK_LIBRARIES}"
      )

    target_link_libraries(${G4LIBTARGET_NAME}-static
      ${G4LIBTARGET_GEANT4_LINK_LIBRARIES_STATIC}
      ${G4LIBTARGET_LINK_LIBRARIES_STATIC})

    # But we can rename the output library to the correct name
    # On WIN32 we *retain* the -static postfix because otherwise
    # we'll conflict with the .lib from the DLL build...
    # We could also install differently...
    if(NOT WIN32)
      set_target_properties(${G4LIBTARGET_NAME}-static
        PROPERTIES OUTPUT_NAME ${G4LIBTARGET_NAME})
    endif()

    set_target_properties(${G4LIBTARGET_NAME}-static
      PROPERTIES CLEAN_DIRECT_OUTPUT 1)

    install(TARGETS ${G4LIBTARGET_NAME}-static
      EXPORT Geant4LibraryDepends
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

    set_property(GLOBAL APPEND
      PROPERTY GEANT4_EXPORTED_TARGETS ${G4LIBTARGET_NAME}-static)
  endif()
endmacro()

#-----------------------------------------------------------------------
# - GEANT4_HEADER_MODULE_TARGET
# Build and install for a header only Geant4 module.
#
macro(geant4_header_module_target)
  CMAKE_PARSE_ARGUMENTS(G4HEADERMOD
    ""
    "COMPONENT"
    ""
    ${ARGN}
    )

  # Only has one component, and we just have to pick out the headers
  include(${G4HEADERMOD_COMPONENT})

  # Header install?
  install(FILES ${${G4MODULENAME}_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    COMPONENT Development)

  # Store the include path of the component so that the build tree
  # config file can pick up all needed header paths
  set_property(GLOBAL APPEND
    PROPERTY GEANT4_BUILDTREE_INCLUDE_DIRS ${${G4MODULENAME}_INCDIR})
endmacro()

#-----------------------------------------------------------------------
# - GEANT4_GRANULAR_LIBRARY_TARGET
# Build and install for a Geant4 module (granular) library
#
macro(geant4_granular_library_target)
  cmake_parse_arguments(G4GRANLIB
    ""
    "COMPONENT"
    ""
    ${ARGN}
    )

  # Granular lib only has one component, but we must pick out
  # the granular dependencies
  include(${G4GRANLIB_COMPONENT})

  # Add the library target, using variables set by the inclusion of
  # the component file
  geant4_library_target(NAME ${G4MODULENAME}
    SOURCES ${${G4MODULENAME}_SOURCES} ${${G4MODULENAME}_HEADERS}
    GEANT4_LINK_LIBRARIES ${${G4MODULENAME}_GRANULAR_DEPENDENCIES}
    LINK_LIBRARIES ${${G4MODULENAME}_LINK_LIBRARIES})

  # Header install?
  install(FILES ${${G4MODULENAME}_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    COMPONENT Development)

  # Store the include path of the component so that the build tree
  # config file can pick up all needed header paths
  set_property(GLOBAL APPEND
    PROPERTY GEANT4_BUILDTREE_INCLUDE_DIRS ${${G4MODULENAME}_INCDIR})
endmacro()

#-----------------------------------------------------------------------
# - GEANT4_GLOBAL_LIBRARY_TARGET
# Build and install of a Geant4 category (global) library
#
macro(geant4_global_library_target)
  cmake_parse_arguments(G4GLOBLIB
    ""
    "NAME"
    "COMPONENTS"
    ${ARGN}
    )

  # We loop over the component sources one at a time,
  # appending properties as we go.
  foreach(_comp ${G4GLOBLIB_COMPONENTS})
    include(${_comp})
    # In case we have a global lib with one component, ensure name gets set
    if(NOT G4GLOBLIB_NAME)
      set(G4GLOBLIB_NAME ${G4MODULENAME})
    endif()

    list(APPEND ${G4GLOBLIB_NAME}_GLOBAL_SOURCES ${${G4MODULENAME}_SOURCES})
    list(APPEND ${G4GLOBLIB_NAME}_GLOBAL_HEADERS ${${G4MODULENAME}_HEADERS})

    list(APPEND ${G4GLOBLIB_NAME}_GLOBAL_DEPENDENCIES
      ${${G4MODULENAME}_GLOBAL_DEPENDENCIES})

    list(APPEND ${G4GLOBLIB_NAME}_LINK_LIBRARIES
      ${${G4MODULENAME}_LINK_LIBRARIES})

    list(APPEND ${G4GLOBLIB_NAME}_BUILDTREE_INCLUDES ${${G4MODULENAME}_INCDIR})
  endforeach()

  # Reset include dirs, setting from module direct dirs only
  # But need to pass this forward to the library target so that
  # it can use this list in target_include_directories...
  #set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "")
  #include_directories(${${G4GLOBLIB_NAME}_BUILDTREE_INCLUDES})

  # Filter out duplicates/self in GLOBAL_DEPENDENCIES and LINK_LIBRARIES
  if(${G4GLOBLIB_NAME}_GLOBAL_DEPENDENCIES)
    list(REMOVE_DUPLICATES ${G4GLOBLIB_NAME}_GLOBAL_DEPENDENCIES)
    list(REMOVE_ITEM ${G4GLOBLIB_NAME}_GLOBAL_DEPENDENCIES ${G4GLOBLIB_NAME})
  endif()
  if(${G4GLOBLIB_NAME}_LINK_LIBRARIES)
    list(REMOVE_DUPLICATES ${G4GLOBLIB_NAME}_LINK_LIBRARIES)
  endif()

  # Now add the library target
  geant4_library_target(NAME ${G4GLOBLIB_NAME}
    SOURCES
      ${${G4GLOBLIB_NAME}_GLOBAL_SOURCES}
      ${${G4GLOBLIB_NAME}_GLOBAL_HEADERS}
    GEANT4_LINK_LIBRARIES
      ${${G4GLOBLIB_NAME}_GLOBAL_DEPENDENCIES}
    LINK_LIBRARIES
      ${${G4GLOBLIB_NAME}_LINK_LIBRARIES})

  # Header install?
  install(FILES ${${G4GLOBLIB_NAME}_GLOBAL_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    COMPONENT Development)

  # Store the include path of the component so that the build tree
  # config file can pick up all needed header paths
  set_property(GLOBAL APPEND
    PROPERTY GEANT4_BUILDTREE_INCLUDE_DIRS ${${G4GLOBLIB_NAME}_BUILDTREE_INCLUDES})
endmacro()

