#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4expat
# Package: Geant4.src.G4externals.G4zlib
#
# CMakeLists.txt for single level library that may be build global or granular
#
# Created on: 03/05/2013
#
# $Date$
# $Revision$
# $Author$
#
#------------------------------------------------------------------------------

include(Geant4MacroLibraryTargets)
if(GEANT4_BUILD_GRANULAR_LIBS)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
else()
    GEANT4_GLOBAL_LIBRARY_TARGET(COMPONENTS sources.cmake)
endif()

# Externals are a special case during migration to new CMake system.
# To ensure compatibility with both old and new systems, need
# to set INCLUDE_DIRECTORIES usage requirement. Does not affect
# build-time inc dirs, simply propagates public includes to
# clients who link to G4zlib target inside Geant4
foreach(__g4zlib_target G4zlib G4zlib-static)
  if(TARGET ${__g4zlib_target})
    target_include_directories(${__g4zlib_target}
      PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
    if(CMAKE_COMPILER_IS_GNUCXX AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "6.99.99"))
       target_compile_options(${__g4zlib_target} PRIVATE "-Wno-implicit-fallthrough")
    endif()
  endif()
endforeach()

