#---Adding all eventgenerator examples subdirectories explicitly 

cmake_minimum_required(VERSION 3.16...3.21)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_CURRENT_SOURCE_DIR}/../../../cmake/Modules ${CMAKE_MODULE_PATH})
find_package(Geant4)
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------

add_subdirectory(particleGun)
add_subdirectory(exgps)
add_subdirectory(userPrimaryGenerator)

 

#----------------------------------------------------------------------------
# decayer6 example requires Pythia6
#
option(PYTHIA6_INTERNAL "Download and compile Pythia6" OFF)
if (NOT PYTHIA6_INTERNAL)
  find_package(Pythia6)
endif()
if(PYTHIA6_FOUND OR PYTHIA6_INTERNAL)
  if (PYTHIA6_INTERNAL)
    message(STATUS "Pythia6 will be downloaded and compiled")
    include(FetchContent)
    FetchContent_Declare(
    p6
    URL      https://pythia.org/download/pythia6/pythia6428-split.tgz
    URL_HASH MD5=5391da11b95aa1d90600bc69da4f84cd
    )
    FetchContent_GetProperties(p6)
    if (NOT p6_POPULATED)
      FetchContent_Populate(p6)
    endif()
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt.pythia6 DESTINATION ${p6_SOURCE_DIR})
    file(RENAME ${p6_SOURCE_DIR}/CMakeLists.txt.pythia6 ${p6_SOURCE_DIR}/CMakeLists.txt)
    add_subdirectory(${p6_SOURCE_DIR} ${p6_BINARY_DIR})
    set(PYTHIA6_LIBRARIES G4Pythia6)
  endif()
  add_subdirectory(pythia/decayer6)
else()
  message(STATUS "G4 Examples: Pythia6 package not found. --> decayer6 example disabled.") 
endif()

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
# py8decayer example requires Pythia8
#
find_package(Pythia8 QUIET)
if (PYTHIA8_FOUND)
  add_subdirectory(pythia/py8decayer)
else()
  message(STATUS "G4 Examples: Pythia8 package not found. --> py8decayer example disabled.") 
endif()
