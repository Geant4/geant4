#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.16...3.21)
project(movies)

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
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(movies movies.cc ${sources} ${headers})
target_link_libraries(movies ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build movies. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.

# All mac files
file(GLOB macs RELATIVE ${PROJECT_SOURCE_DIR} *mac)
foreach(_mac ${macs})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_mac}
    ${PROJECT_BINARY_DIR}/${_mac}
    COPYONLY
)
endforeach()

# All the example viewfiles
file(GLOB viewfiles RELATIVE ${PROJECT_SOURCE_DIR} viewfiles/*/*)
foreach(_file ${viewfiles})
#  message(STATUS "file=${_file}")
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_file}
    ${PROJECT_BINARY_DIR}/${_file}
    COPYONLY
)
endforeach()

# The old visTutor files
file(GLOB visTutorfiles RELATIVE ${PROJECT_SOURCE_DIR} visTutor/*)
foreach(_file ${visTutorfiles})
#  message(STATUS "file=${_file}")
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_file}
    ${PROJECT_BINARY_DIR}/${_file}
    COPYONLY
)
endforeach()

set(other_files movies.in movies.out movies.err)
foreach(_file ${other_files})
#  message(STATUS "file=${_file}")
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_file}
    ${PROJECT_BINARY_DIR}/${_file}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS movies DESTINATION bin)
