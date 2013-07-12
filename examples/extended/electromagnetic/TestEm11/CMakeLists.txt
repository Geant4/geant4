#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(TestEm11)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Find HBOOK (optional package)
#
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../common/analysis/cmake ${CMAKE_MODULE_PATH})
find_package(HBOOK QUIET)
if(NOT HBOOK_FOUND)
  message(STATUS "G4 Examples: HBOOK package not found. --> TestEm11: g4tools/hbook analysis disabled")  
else()
  message(STATUS "G4 Examples: HBOOK package found. --> TestEm11: g4tools/hbook analysis enabled")  
  add_definitions(-DG4_USE_HBOOK)  
endif()

#----------------------------------------------------------------------------
# Add g4tools/hbook sources
#

if (HBOOK_FOUND)
  include_directories(
       ${PROJECT_SOURCE_DIR}/../../common/analysis/include)

  file(GLOB add_sources ${PROJECT_SOURCE_DIR}/../../common/analysis/src/*.cc)
  list(APPEND sources ${add_sources})

  add_custom_command(
       OUTPUT  ${PROJECT_BINARY_DIR}/close.o
       COMMAND gfortran
       ARGS  -c ${Geant4_INCLUDE_DIR}/tools/hbook/close.f )

  add_custom_command(
       OUTPUT  ${PROJECT_BINARY_DIR}/setpawc.o
       COMMAND gfortran
       ARGS  -c ${Geant4_INCLUDE_DIR}/tools/hbook/setpawc.f )

  add_custom_command(
       OUTPUT  ${PROJECT_BINARY_DIR}/setntuc.o
       COMMAND gfortran
       ARGS  -c ${Geant4_INCLUDE_DIR}/tools/hbook/setntuc.f )
  set(TOOLS_FORTRAN_OBJECTS close.o setpawc.o setntuc.o)       
endif()       

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(TestEm11 TestEm11.cc ${sources} ${headers} ${TOOLS_FORTRAN_OBJECTS})
target_link_libraries(TestEm11 ${Geant4_LIBRARIES} ${HBOOK_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build TestEm11. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(TestEm11_SCRIPTS
    multiLayers.mac range.mac run01.mac run02.mac sandia.mac TestEm11.in TestEm11.out vis.mac water.mac
  )
  
foreach(_script ${TestEm11_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS TestEm11 DESTINATION bin)

