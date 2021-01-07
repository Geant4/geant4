#---Adding all eventgenerator examples subdirectories explicitly 

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

find_package(Geant4)
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------

add_subdirectory(particleGun)
add_subdirectory(exgps)
add_subdirectory(userPrimaryGenerator)

#----------------------------------------------------------------------------
# HepMC examples require HepMC
#
find_package(HepMC QUIET)
if(HEPMC_FOUND)
  add_subdirectory(HepMC)
else()
  message(STATUS "G4 Examples: HepMC package not found. --> HepMC examples disabled.") 
endif()   

#----------------------------------------------------------------------------
# decayer6 example requires Pythia6
#
find_package(Pythia6 QUIET)
if(PYTHIA6_FOUND)
  add_subdirectory(pythia/decayer6)
else()
  message(STATUS "G4 Examples: Pythia6 package not found. --> decayer6 example disabled.") 
endif()   
