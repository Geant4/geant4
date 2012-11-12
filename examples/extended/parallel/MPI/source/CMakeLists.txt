# - CmakeLists.txt for MPI interface

#------------------------------------------------------------------------------
project(G4MPI)
#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

# check MPI package...
find_package(MPI REQUIRED)
# modify these var. if needed
set(CMAKE_CXX_COMPILER mpicxx)
#set(CMAKE_CXX_INCLUDE_PATH )

#------------------------------------------------------------------------------
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})

# sources
add_library(G4mpi SHARED
  src/G4MPIbatch.cc
  src/G4MPImanager.cc
  src/G4MPImessenger.cc
  src/G4MPIrandomSeedGenerator.cc
  src/G4MPIsession.cc
  src/G4MPIstatus.cc
  src/G4UImpish.cc
  src/G4VMPIseedGenerator.cc
  src/G4VMPIsession.cc
)

# headers
set(HEADERS
  include/G4MPIbatch.hh
  include/G4MPImanager.hh
  include/G4MPImessenger.hh
  include/G4MPIrandomSeedGenerator.hh
  include/G4MPIsession.hh
  include/G4MPIstatus.hh
  include/G4UImpish.hh
  include/G4VMPIseedGenerator.hh
  include/G4VMPIsession.hh
)

# install
INSTALL(TARGETS G4mpi LIBRARY DESTINATION ${Geant4_DIR}/..)
INSTALL(FILES ${HEADERS} DESTINATION ${Geant4_INCLUDE_DIR})

