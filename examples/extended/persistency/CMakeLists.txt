#---Adding all persistency examples subdirectories explicitly 

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

find_package(Geant4)
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# P01, P02 examples require ROOT, GCCXML, BUILD_SHARED_LIBS
#
find_package(ROOT QUIET)
if(ROOT_FOUND AND Geant4_shared_FOUND)
  #if(ROOT_FOUND AND GCCXML)
  STRING(REGEX MATCH "6.*" VERSION6MATCH ${ROOT_VERSION})
  if(NOT VERSION6MATCH) 
    message("-- G4 Examples : This example requires root6, disabled")
    return()
  else()
    add_subdirectory(P01)
    add_subdirectory(P02)
  endif()
else()
  if(NOT ROOT_FOUND)
    message(STATUS "G4 Examples: ROOT package not found. ") 
  endif()   
  if(NOT Geant4_shared_FOUND)
    message(STATUS "G4 Examples: Shared Libraires not available.")
  endif()   
  message(STATUS "P01, P02 disabled")
endif()

#----------------------------------------------------------------------------
# gdml examples require Geant4 build with gdml
#
if(Geant4_gdml_FOUND)
  add_subdirectory(gdml)
else()  
  message(STATUS "G4 Examples: GDML not installed --> gdml examples disabled")  
endif()
