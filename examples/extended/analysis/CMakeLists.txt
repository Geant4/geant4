#---Adding all analysis examples subdirectories explicitly 

cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

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
