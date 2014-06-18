# - CmakeLists.txt for MPI interface

#------------------------------------------------------------------------------
project(G4MPI)
#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

# check MPI package...
find_package(MPI REQUIRED)

#------------------------------------------------------------------------------
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

#if(NOT CMAKE_BUILD_TYPE)
#  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
#endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})

# sources
add_library(G4mpi
  src/G4MPIbatch.cc
  src/G4MPImanager.cc
  src/G4MPImessenger.cc
  src/G4MPIrandomSeedGenerator.cc
  src/G4MPIsession.cc
  src/G4MPIstatus.cc
  src/G4UImpish.cc
  src/G4VMPIseedGenerator.cc
  src/G4VMPIsession.cc
  src/G4MPIScorerMerger.cc
  src/G4MPIRunMerger.cc
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
  include/G4MPIScorerMerger.hh
  include/G4MPIRunMerger.hh
)

# install
INSTALL(TARGETS G4mpi DESTINATION lib)
INSTALL(FILES ${HEADERS} DESTINATION include)
INSTALL(FILES ${PROJECT_BINARY_DIR}/G4mpiConfig.cmake DESTINATION lib/G4mpi)

#Generate G4mpiConfig.cmake file
#Get full path of library in installation directory
GET_PROPERTY( G4MPIFULLPATH TARGET G4mpi PROPERTY LOCATION)
GET_FILENAME_COMPONENT(G4MPIFULLPATH ${G4MPIFULLPATH} NAME)
SET( G4MPIFULLPATH "${CMAKE_INSTALL_PREFIX}/lib/${G4MPIFULLPATH}")
SET(G4MPILIB G4mpi)
CONFIGURE_FILE(G4mpiConfig.cmake.in
			  ${PROJECT_BINARY_DIR}/G4mpiConfig.cmake
			  @ONLY)

