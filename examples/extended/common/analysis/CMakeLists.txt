#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(analysis)

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find HBOOK (required package)
#
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
find_package(HBOOK QUIET)
if(NOT HBOOK_FOUND)
  message(STATUS "G4 Examples: HBOOK package not found. --> g4tools/hbook analysis disabled")  
  return()
else()
  message(STATUS "G4 Examples: HBOOK package found. --> g4tools/hbook analysis enabled")  
  add_definitions(-DG4_USE_HBOOK)  
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Build fortran sources from g4tools/hbook
#

if(HBOOK_FOUND)
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
add_executable(testAnalysis testAnalysis.cc ${sources} ${headers} ${TOOLS_FORTRAN_OBJECTS})
target_link_libraries(testAnalysis ${Geant4_LIBRARIES} ${HBOOK_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build analysis. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(analysis_SCRIPTS
   
  )

foreach(_script ${analysis_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(analysis DEPENDS testAnalysis)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS testAnalysis DESTINATION bin)

