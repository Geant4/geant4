#Workaround: rely on executables build by cmake since it does not compile
#correctly because FindG4mpi.cmake is not found
find_package(MPI QUIET)
if(NOT MPI_CXX_FOUND)
  message(STATUS "G4 Examples: MPI not found --> mpi based tests disabled")
  return()
endif()
message(STATUS "G4 Examples: mpi examples will use mpi launcher: " ${MPIEXEC})

# Base output dir for MPI examples:
set(G4MPI_CTESTS_BASE_OUTPUT_DIR "${PROJECT_BINARY_DIR}/examples/extended/parallel/MPI")

# Set G4mpi_DIR fo later pass-down to examples
set(G4mpi_DIR "${G4MPI_CTESTS_BASE_OUTPUT_DIR}/G4mpi")

# - Build/Test G4mpi
geant4_add_test(mpi-libg4mpi
  COMMAND ${CMAKE_COMMAND} -E echo "G4mpi build complete"
  SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/source
  BINARY_DIR ${G4mpi_DIR}
  PROJECT libG4mpi
  BUILD G4mpi
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  LABELS MPI
  )

# - Build/Test exMPI01
# variable to simplify paths
set(G4MPI_EX01_BINDIR ${G4MPI_CTESTS_BASE_OUTPUT_DIR}/exMPI01)

geant4_add_test(mpi-ex01-sequential
  SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/examples/exMPI01
  BINARY_DIR ${G4MPI_EX01_BINDIR}
  PROJECT exMPI01
  BUILD exMPI01
  COMMAND ${G4MPI_EX01_BINDIR}/exMPI01 run.mac
  WORKING_DIRECTORY ${G4MPI_EX01_BINDIR}
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  DEPENDS mpi-libg4mpi
  LABELS MPI
  )
# Needs G4mpi
set_property(TEST mpi-ex01-sequential APPEND PROPERTY ENVIRONMENT G4mpi_DIR=${G4mpi_DIR})

geant4_add_test(mpi-ex01
  COMMAND ${MPIEXEC} -n 2 ${G4MPI_EX01_BINDIR}/exMPI01 run.mac
  WORKING_DIRECTORY ${G4MPI_EX01_BINDIR}
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  DEPENDS mpi-ex01-sequential
  LABELS MPI
  )

# - Build/Test exMPI02
# Needs ROOT
find_package(ROOT QUIET)
if(ROOT_FOUND)
  set(G4MPI_EX02_BINDIR ${G4MPI_CTESTS_BASE_OUTPUT_DIR}/exMPI02)

  geant4_add_test(mpi-ex02-sequential
    SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/examples/exMPI02
    BINARY_DIR ${G4MPI_EX02_BINDIR}
    PROJECT exMPI02
    BUILD exMPI02
    COMMAND ${G4MPI_EX02_BINDIR}/exMPI02 run.mac
    WORKING_DIRECTORY ${G4MPI_EX02_BINDIR}
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS mpi-libg4mpi
    LABELS MPI
    )
  # Needs G4mpi
  set_property(TEST mpi-ex02-sequential APPEND PROPERTY ENVIRONMENT G4mpi_DIR=${G4mpi_DIR})

  geant4_add_test(mpi-ex02
    COMMAND ${MPIEXEC} -n 2 ${G4MPI_EX02_BINDIR}/exMPI02 run.mac
    WORKING_DIRECTORY ${G4MPI_EX02_BINDIR}
		ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
		DEPENDS mpi-ex02-sequential
		LABELS MPI
    )
else()
    message(STATUS "G4 Examples: MPI example exMPI02 disabled: ROOT not found")
endif()

# - Build/Test exMPI03
# variable to simplify paths
set(G4MPI_EX03_BINDIR ${G4MPI_CTESTS_BASE_OUTPUT_DIR}/exMPI03)

geant4_add_test(mpi-ex03-sequential
  SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/examples/exMPI03
  BINARY_DIR ${G4MPI_EX03_BINDIR}
  PROJECT exMPI03
  BUILD exMPI03
  COMMAND ${G4MPI_EX03_BINDIR}/exMPI03 run.mac
  WORKING_DIRECTORY ${G4MPI_EX03_BINDIR}
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  DEPENDS mpi-libg4mpi
  LABELS MPI
  )
# Needs G4mpi
set_property(TEST mpi-ex03-sequential APPEND PROPERTY ENVIRONMENT G4mpi_DIR=${G4mpi_DIR})

geant4_add_test(mpi-ex03
  COMMAND ${MPIEXEC} -n 2 ${G4MPI_EX03_BINDIR}/exMPI03 run.mac
  WORKING_DIRECTORY ${G4MPI_EX03_BINDIR}
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  DEPENDS mpi-ex03-sequential
  LABELS MPI
  )

# - Build/Test exMPI04
# variable to simplify paths
set(G4MPI_EX04_BINDIR ${G4MPI_CTESTS_BASE_OUTPUT_DIR}/exMPI04)

geant4_add_test(mpi-ex04-sequential
  SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/examples/exMPI04
  BINARY_DIR ${G4MPI_EX04_BINDIR}
  PROJECT exMPI04
  BUILD exMPI04
  COMMAND ${G4MPI_EX04_BINDIR}/exMPI04 run.mac
  WORKING_DIRECTORY ${G4MPI_EX04_BINDIR}
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  DEPENDS mpi-libg4mpi
  LABELS MPI
  )
# Needs G4mpi
set_property(TEST mpi-ex04-sequential APPEND PROPERTY ENVIRONMENT G4mpi_DIR=${G4mpi_DIR})

geant4_add_test(mpi-ex04
  COMMAND ${MPIEXEC} -n 3 ${G4MPI_EX04_BINDIR}/exMPI04 run.mac
  WORKING_DIRECTORY ${G4MPI_EX04_BINDIR}
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  DEPENDS mpi-ex04-sequential
  LABELS MPI
  )
