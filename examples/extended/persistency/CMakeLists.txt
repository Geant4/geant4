#---Adding all persistency examples subdirectories explicitly 

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

add_subdirectory(P01)
add_subdirectory(P02)
add_subdirectory(P03)

#----------------------------------------------------------------------------
# gdml examles require Geant4 build with gdml
#
find_package(Geant4)
if(Geant4_gdml_FOUND)
  add_subdirectory(gdml)
else()  
  message(STATUS "G4 Examples: GDML not installed --> gdml examples disabled")  
endif()
