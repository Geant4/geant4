#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(wholeNuclearDNA)

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
# Add the executable, and link it to the Geant4 libraries
#

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cc ${sources} ${headers})
target_link_libraries(${PROJECT_NAME} ${Geant4_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build microdosimetry. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(dnageometry_SCRIPTS
    wholenucleardna.in 
    vis.mac 
    chromo1.dat 
    chromo2.dat 
    chromo3.dat 
    chromo4.dat 
    chromo5.dat 
    chromo6.dat 
    chromo7.dat 
    chromo8.dat 
    chromo9.dat 
    chromo10.dat 
    chromo11.dat 
    chromo12.dat 
    chromo13.dat 
    chromo14.dat 
    chromo15.dat 
    chromo16.dat 
    chromo17.dat 
    chromo18.dat 
    chromo19.dat 
    chromo20.dat 
    chromo21.dat 
    chromo22.dat 
    chromoX.dat 
    chromoY.dat 
  )

file(GLOB SCRIPTS ${PROJECT_SOURCE_DIR}/*.mac 
                  ${PROJECT_SOURCE_DIR}/*.in
                  ${PROJECT_SOURCE_DIR}/*.dat
                  ${PROJECT_SOURCE_DIR}/*.C 
                    )

foreach(_script ${SCRIPTS})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${PROJECT_NAME} DESTINATION bin)

