#------------------------------------------------------------------------------
# Define library
# We only build as a static library because we only want to use this in the 
# example applications. It also avoids issues with dylib/DLL export/import
# and paths on some platforms.
add_library(G4mpi STATIC
  include/G4MPIbatch.hh
  include/G4MPIextraWorker.hh
  include/G4MPImanager.hh
  include/G4MPImessenger.hh
  include/G4MPIrandomSeedGenerator.hh
  include/G4MPIsession.hh
  include/G4MPIstatus.hh
  include/G4UImpish.hh
  include/G4VMPIextraWorker.hh
  include/G4VMPIseedGenerator.hh
  include/G4VMPIsession.hh
  include/G4MPIscorerMerger.hh
  include/G4MPIrunMerger.hh
  include/G4MPIhistoMerger.hh
  include/G4MPIntupleMerger.hh
  include/G4VUserMPIrunMerger.hh
  include/G4MPIutils.hh
  src/G4MPIbatch.cc
  src/G4MPIextraWorker.cc
  src/G4MPImanager.cc
  src/G4MPImessenger.cc
  src/G4MPIrandomSeedGenerator.cc
  src/G4MPIsession.cc
  src/G4MPIstatus.cc
  src/G4UImpish.cc
  src/G4VMPIseedGenerator.cc
  src/G4VMPIsession.cc
  src/G4MPIscorerMerger.cc
  src/G4MPIhistoMerger.cc
  src/G4MPIntupleMerger.cc
  src/G4VUserMPIrunMerger.cc
  src/G4MPIutils.cc
  analysis/include/G4RootMpiAnalysisManager.hh
  analysis/include/G4RootMpiNtupleFileManager.hh
  analysis/include/G4RootMpiNtupleManager.hh
  analysis/include/G4RootMpiPNtupleDescription.hh
  analysis/include/G4RootMpiPNtupleManager.hh
  analysis/src/G4RootMpiAnalysisManager.cc
  analysis/src/G4RootMpiNtupleFileManager.cc
  analysis/src/G4RootMpiNtupleManager.cc
  analysis/src/G4RootMpiPNtupleManager.cc
)
target_include_directories(G4mpi 
  PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/analysis/include>
)
target_link_libraries(G4mpi
  PUBLIC
    MPI::MPI_CXX 
    Geant4::G4analysis
    Geant4::G4intercoms
    Geant4::G4interfaces
    Geant4::G4global
    Geant4::G4run
    Geant4::G4tools
)

#------------------------------------------------------------------------------
# Create a namespaced alias target for the library for forward compatibility
# in case the library is moved or reorganised in the future
add_library(G4mpi::G4mpi ALIAS G4mpi)
