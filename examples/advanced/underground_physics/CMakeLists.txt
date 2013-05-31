#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(underground_physics)

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

#Check here for DMXENV_GPS_USE
#Use GPS if the environment variable is set, and the user has *not* supplied the CMake 
#command line argument. Otherwise, the argument of the command line has the priority.
#
if("$ENV{DMXENV_GPS_USE}" AND NOT DEFINED DMXENV_GPS_USE)
  set(DMXENV_GPS_USE_DEFAULT ON)
else()
  set(DMXENV_GPS_USE_DEFAULT ${DMXENV_GPS_USE})
endif()

option(DMXENV_GPS_USE "Build DMX with support for GPS" ${DMXENV_GPS_USE_DEFAULT})
mark_as_advanced(DMXENV_GPS_USE)

if(DMXENV_GPS_USE)
  add_definitions(-DDMXENV_GPS_USE)
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
add_executable(DMX DMX.cc ${sources} ${headers})
target_link_libraries(DMX ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build underground_physics. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(underground_physics_SCRIPTS
    alpha_1000.mac alpha_low.mac alpha.mac ambe_spectrum.mac gamma_1000.mac gamma_2000.mac gamma_low.mac gamma.mac initInter.mac neutron.mac stt.mac stt.out stt_pmt_ref.out stt_ref.out wired.mac
  )

foreach(_script ${underground_physics_SCRIPTS})
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
add_custom_target(underground_physics DEPENDS DMX)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS DMX DESTINATION bin)

