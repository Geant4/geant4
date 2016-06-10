#---Adding all analysis examples subdirectories explicitly 

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

find_package(Geant4)
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------

add_subdirectory(AnaEx01)
add_subdirectory(B1Con)

#----------------------------------------------------------------------------
# AnaEx02 example requires ROOT
#
find_package(ROOT QUIET)
if(ROOT_FOUND)
  add_subdirectory(AnaEx02)
else()
  message(STATUS "G4 Examples: ROOT package not found. AnaEx02 disabled.") 
endif()   

#----------------------------------------------------------------------------
# AnaEx03 example requires AIDA
#
find_package(AIDA QUIET)
if(AIDA_FOUND)
  add_subdirectory(AnaEx03)
else()
  message(STATUS "G4 Examples: AIDA package not found. AnaEx03 disabled.") 
endif()   
