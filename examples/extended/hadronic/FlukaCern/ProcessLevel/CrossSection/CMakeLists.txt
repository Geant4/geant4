#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(HadronNucleusXS)

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
# Check whether FLUKA should be used or not
set(G4_USE_FLUKA OFF CACHE BOOL "Using FLUKA")
if(G4_USE_FLUKA)
  message(STATUS "G4_USE_FLUKA=ON : Using FLUKA interface for building ${PROJECT_SOURCE_DIR}")
  add_definitions(-DG4_USE_FLUKA)
  find_package(FLUKAInterface REQUIRED)
  if(FLUKAInterface_FOUND)
    message(STATUS "FLUKA cmake module was found : ${CMAKE_MODULE_PATH}")
  else()
    message(FATAL_ERROR "FLUKA cmake module was NOT found! Please add one.")
  endif()
else()
  message(STATUS "G4_USE_FLUKA=OFF : NOT using FLUKA interface for building ${PROJECT_SOURCE_DIR}. \n \
  If ever you want to use the FLUKA interface, please repeat cmake command with -DG4_USE_FLUKA=1")
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${PROJECT_SOURCE_DIR}/../../utils/include
                    ${FLUKAInterface_INCLUDE_DIR}
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/../../utils/src/*.cc ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/../../utils/include/*.hh ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(HadronNucleusXS HadronNucleusXS.cc ${sources} ${headers})
target_link_libraries(HadronNucleusXS ${FLUKAInterface_LIBRARIES} ${Geant4_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build HadronNucleusXS. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(HadronNucleusXS_SCRIPTS
    all_XS.in
  )

foreach(_script ${HadronNucleusXS_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS HadronNucleusXS DESTINATION bin)
