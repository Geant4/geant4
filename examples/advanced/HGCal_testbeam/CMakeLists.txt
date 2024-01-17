#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)

project(HGCal_testbeam)

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
# Find ROOT (required package only for reading input particles from file)
#
find_package(ROOT QUIET NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
if (ROOT_FOUND)
  message(STATUS "Found ROOT: ${ROOT_DIR}. Particles read from file can be used as generator.")
  add_definitions(-DWITHROOT)
  include_directories(${ROOT_INCLUDE_DIRS})
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
add_executable(HGCal_testbeam hgcal_testbeam.cc ${sources} ${headers})
target_link_libraries(HGCal_testbeam ${Geant4_LIBRARIES})
if (ROOT_FOUND)
  target_link_libraries(HGCal_testbeam ${ROOT_LIBRARIES})
  set_target_properties(HGCal_testbeam PROPERTIES CXX_STANDARD 17)
endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build HGCal)testbeam. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(HGCAL_SCRIPTS
  init_vis.mac
  run.mac
  vis.mac
  )

foreach(_script ${HGCAL_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

if (ROOT_FOUND)
  configure_file(
    ${PROJECT_SOURCE_DIR}/readFromFile.mac
    ${PROJECT_BINARY_DIR}/readFromFile.mac
    COPYONLY
    )
endif()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS HGCal_testbeam DESTINATION bin)


