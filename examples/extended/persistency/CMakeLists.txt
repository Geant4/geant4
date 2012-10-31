#---Adding all persistency examples subdirectories explicitly 

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

if(GEANT4_USE_GDML)
  add_subdirectory(gdml)
endif()  
add_subdirectory(P01)
add_subdirectory(P02)
add_subdirectory(P03)
