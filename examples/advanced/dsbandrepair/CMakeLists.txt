cmake_minimum_required(VERSION 3.16...3.27)
project(dsbandrepair)

find_package(Geant4 REQUIRED)

option(DOWNLOAD_GEOMETRY "Download geometry files" TRUE)
option(USE_MPI "Using MPI" FALSE)
if (USE_MPI)
  find_package(G4mpi REQUIRED)
endif()
#----------------------------------------------------------------------------
if (DOWNLOAD_GEOMETRY) 
  include(ExternalProject)
  ExternalProject_Add(dnafabric_geometries
  SOURCE_DIR ${PROJECT_BINARY_DIR}/dnafabric_geometries
  URL https://cern.ch/geant4-data/datasets/examples/advanced/dsbandrepair/0/dnafabric_geometries.tar.xz
  URL_HASH SHA256=7e77ec0dd4291599768a4c95b95f2456a1477b32f8141c97ddb78f982d828649
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  )
else()
  message("---> Option for downloading geometry files is off. Make sure you already have them, or download manually!")
endif()
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(${PROJECT_NAME} dsbandrepair.cc ${sources} ${headers})
target_link_libraries(${PROJECT_NAME} ${Geant4_LIBRARIES} ${G4mpi_LIBRARIES})
#----------------------------------------------------------------------------
if(DOWNLOAD_GEOMETRY)
  add_dependencies(${PROJECT_NAME} dnafabric_geometries)
endif()
#----------------------------------------------------------------------------
if (USE_MPI)
  message(STATUS "dsbandrepair will run with MPI")
  target_compile_definitions(${PROJECT_NAME} PRIVATE USE_MPI)
endif()
#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# 
target_include_directories(${PROJECT_NAME} PUBLIC
                    ${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR}
                    ${G4mpi_INCLUDE_DIR}
)
#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build dsbandrepair. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(dsbandrepair_SCRIPTS
  ${PROJECT_SOURCE_DIR}/macros/dsbandrepair.in
  ${PROJECT_SOURCE_DIR}/macros/chem.in
  ${PROJECT_SOURCE_DIR}/macros/endophys.in
  ${PROJECT_SOURCE_DIR}/macros/fibroblast.in
  ${PROJECT_SOURCE_DIR}/macros/yeastphys.in
  ${PROJECT_SOURCE_DIR}/macros/analysis.in
  )

foreach(_script ${dsbandrepair_SCRIPTS})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
    )
endforeach()
