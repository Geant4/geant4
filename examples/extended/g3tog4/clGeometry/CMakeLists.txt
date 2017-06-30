#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(clGeometry)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED g3tog4 ui_all vis_all)
else()
  find_package(Geant4 REQUIRED g3tog4)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/common/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add common subdirectory
#
set (common_prefix clGeometry_)
add_subdirectory(common) 

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(clGeometry clGeometry.cc ${sources} ${headers})
target_link_libraries(clGeometry ${common_prefix}common ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build clGeometry. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(clGeometry_SCRIPTS
    clGeometry.in clGeometry.out init.mac init_vis.mac vis.mac
  )

foreach(_script ${clGeometry_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Copy data
#
file(GLOB clGeometry_DATA 
     RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/data/*.dat)
foreach(_data_file ${clGeometry_DATA})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_data_file}
    ${PROJECT_BINARY_DIR}/${_data_file}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS clGeometry DESTINATION bin)

