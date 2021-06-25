#---Adding all parameterisations examples subdirectories explicitly 

cmake_minimum_required(VERSION 3.12...3.20)

find_package(Geant4)
include(${Geant4_USE_FILE})

add_subdirectory(gflash)
add_subdirectory(Par01)
add_subdirectory(Par03)

#----------------------------------------------------------------------------
# Par02 example requires Geant4 build with gdml
#
if(Geant4_gdml_FOUND)
  add_subdirectory(Par02)
else()  
  message(STATUS "G4 Examples: GDML not installed --> Par02 example disabled")  
endif()
