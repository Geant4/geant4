#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(${common_prefix}common)

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
find_package(Geant4 REQUIRED)

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
# Add the library
#
add_library(${common_prefix}common ${sources} ${headers})
target_link_libraries(${common_prefix}common ${Geant4_LIBRARIES} )

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries and common library
#
add_executable(${common_prefix}testCommon testCommon.cc )
target_link_libraries(${common_prefix}testCommon ${common_prefix}common ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build detectorConstruction. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(detectorConstruction_SCRIPTS
   
  )

foreach(_script ${detectorConstruction_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
# and the library to 'lib' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${common_prefix}testCommon DESTINATION bin)
install(TARGETS ${common_prefix}common DESTINATION lib)
