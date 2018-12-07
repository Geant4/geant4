#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(hadrontherapy)

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
# Find ROOT (required package)
#
# find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${PROJECT_SOURCE_DIR}/shared/include
                    ${Geant4_INCLUDE_DIR})
		    # ${ROOT_INCLUDE_DIRS})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc
                  ${PROJECT_SOURCE_DIR}/shared/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh
                  ${PROJECT_SOURCE_DIR}/shared/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(hadrontherapy hadrontherapy.cc ${sources} ${headers})
target_link_libraries(hadrontherapy ${Geant4_LIBRARIES})
# ${ROOT_LIBRARIES} )
file(GLOB_RECURSE dataFiles RELATIVE ${PROJECT_SOURCE_DIR} data/rbe/*.*)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build hadrontherapy. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
file(GLOB macros RELATIVE ${PROJECT_SOURCE_DIR} macro/*.mac  macro/*.loop)
file(GLOB modulators RELATIVE ${PROJECT_SOURCE_DIR} Modulators/*.txt)
file(GLOB fields RELATIVE ${PROJECT_SOURCE_DIR} field/*.TABLE)
file(GLOB_RECURSE data RELATIVE ${PROJECT_SOURCE_DIR} data/*)

set(other_files
  clean.sh
  batch.mac
)

foreach(_file ${macros} ${modulators} ${fields} ${other_files} ${data})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_file}
    ${PROJECT_BINARY_DIR}/${_file}
    COPYONLY
    )
endforeach()


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS hadrontherapy DESTINATION bin)

