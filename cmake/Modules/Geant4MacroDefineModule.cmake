# This file defines the following macro for Geant4 developers needing to
# define the sources, headers and library dependencies for a standard 
# Geant4 granular library module:
#
# GEANT4_DEFINE_MODULE      - define a standard Geant4 Granular Library Module
#
# A Geant4 Module is defined as a directory containing subdirectories
#
#   include - holds all header files for the module
#   src     - holds all source files for the module
#
# GEANT4_DEFINE_MODULE will take the name of the module, a list of header files,
# a list of source files and dependencies:
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
# headers and sources, defining the variables
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

include(CMakeMacroParseArguments)

MACRO(GEANT4_DEFINE_MODULE)
    PARSE_ARGUMENTS(G4DEFMOD
        "NAME;HEADERS;SOURCES;GRANULAR_DEPENDENCIES;GLOBAL_DEPENDENCIES;LINK_LIBRARIES"    
        ""
        ${ARGN}
    )

    set(G4MODULENAME ${G4DEFMOD_NAME})

    get_filename_component(${G4MODULENAME}_BASEDIR ${CMAKE_CURRENT_LIST_FILE} PATH)
    set(${G4MODULENAME}_SRCDIR ${${G4MODULENAME}_BASEDIR}/src)
    set(${G4MODULENAME}_INCDIR ${${G4MODULENAME}_BASEDIR}/include)

    foreach(_HDR ${G4DEFMOD_HEADERS})
        list(APPEND ${G4MODULENAME}_HEADERS ${${G4MODULENAME}_INCDIR}/${_HDR})
    endforeach()


    foreach(_SRC ${G4DEFMOD_SOURCES})
        list(APPEND ${G4MODULENAME}_SOURCES ${${G4MODULENAME}_SRCDIR}/${_SRC})
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
ENDMACRO()

