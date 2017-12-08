# - CmakeLists.txt for building an application

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(exMPI02)

#------------------------------------------------------------------------------
# check MPI package...
find_package(MPI REQUIRED)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

find_package(G4mpi REQUIRED)
#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
  message(STATUS "G4 Examples: ROOT package not found. --> exMPI02 example disabled")
  return()
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIRS}
		    ${G4mpi_INCLUDE_DIR})


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exMPI02 exMPI02.cc ${sources} ${headers})
target_link_libraries(exMPI02 ${G4mpi_LIBRARIES} ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build exMPI02. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(exMPI02_SCRIPTS
    run.mac
    v.mac
  )

foreach(_script ${exMPI02_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS exMPI02 RUNTIME DESTINATION bin)
