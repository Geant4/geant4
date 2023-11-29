#---Adding all analysis examples subdirectories explicitly

cmake_minimum_required(VERSION 3.16...3.21)

find_package(Geant4)
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------

add_subdirectory(AnaEx01)
add_subdirectory(AnaEx03)
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
