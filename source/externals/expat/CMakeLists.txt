#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4expat
# Package: Geant4.src.G4externals.G4expat
#
# CMakeLists.txt for single level library that may be build global or granular
#
# Created on: 02/06/2011
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

# - Ensure clients of the static expat library use the XML_STATIC
# preprocessor def
if(TARGET G4expat-static)
  target_compile_definitions(G4expat-static PUBLIC XML_STATIC)
endif()

# Externals are a special case during migration to new CMake system.
# To ensure compatibility with both old and new systems, need
# to set INCLUDE_DIRECTORIES usage requirement. Does not affect
# build-time inc dirs, simply propagates public includes to
# clients who link to G4expat target inside Geant4
foreach(__g4expat_target G4expat G4expat-static)
  if(TARGET ${__g4expat_target})
    target_include_directories(${__g4expat_target}
      PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
      )
  endif()
endforeach()

