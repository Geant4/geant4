#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(dnadamage1)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available Vis drivers by default
# You can set WITH_GEANT4_VIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_VIS "Build example with Geant4 Vis drivers" ON)
if(WITH_GEANT4_VIS)
  find_package(Geant4 REQUIRED vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Dowload geometry data file

set(GEOMETRY_NEEDS_DOWNLOAD TRUE)
set(GEOMETRY_FILE_NAME "VoxelStraight.fab2g4dna")
set(GEOMETRY_LOCAL_FILENAME "${PROJECT_BINARY_DIR}/${GEOMETRY_FILE_NAME}")
set(GEOMETRY_DATASETS_URL 
"https://cern.ch/geant4-data/examples/dna/dnadamage1/${GEOMETRY_FILE_NAME}")
set(HASH_MD5 "3e28151dc4c4647af3ae37d0385fc443")

if(EXISTS "${GEOMETRY_FILE_NAME}")
  set(GEOMETRY_NEEDS_DOWNLOAD FALSE)
endif()


if(GEOMETRY_NEEDS_DOWNLOAD)
  message(STATUS "Geometry-data: attempting download: ${GEOMETRY_DATASETS_URL}")
  file(DOWNLOAD "${GEOMETRY_DATASETS_URL}" "${GEOMETRY_LOCAL_FILENAME}"
    INACTIVITY_TIMEOUT 50
    TIMEOUT 50
    STATUS DownloadStatus
    )

  list(GET DownloadStatus 0 DownloadReturnStatus)
  if(DownloadReturnStatus)
    message(STATUS "Geometry-data: download FAILED: ${DownloadReturnStatus}, 
    This example needs internet for the geometry data file,
    even configuring done and complied. 
    Please, check your connection.
    ")
  else()
    message(STATUS "Geometry-data: download OK")
  endif()
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#

include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(dnadamage1 dnadamage1.cc ${sources} ${headers})
target_link_libraries(dnadamage1 ${Geant4_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build dnadamage1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

file(GLOB MAC_FILES ${PROJECT_SOURCE_DIR}/*.mac 
     ${PROJECT_SOURCE_DIR}/*.in
     ${PROJECT_SOURCE_DIR}/*.C
	${PROJECT_SOURCE_DIR}/*.fab2g4dna)

foreach(_script ${MAC_FILES})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS dnadamage1 DESTINATION bin)


