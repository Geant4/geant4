#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(radioprotection)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available Vis drivers by default
# You can set WITH_GEANT4_VIS to OFF via the command line or ccmake/cmake-gui
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
# Setup of ROOT ANALYSIS : optional.
option(WITH_ANALYSIS_USE "Build example with analysis objects" ON)
if(WITH_ANALYSIS_USE)
add_definitions(-DANALYSIS_USE)
else()
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                     ${Geant4_INCLUDE_DIR} )
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(radioprotection radioprotection.cc ${sources} ${headers})
target_link_libraries(radioprotection ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build basic_radioprotection. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(radioprotection_SCRIPTS
   vis.mac physics.mac primary.mac run.mac geometry.mac output.mac  MergeFiles LoadPlotNtuple.C ProcessMicro.C 1_plot_distributions.py 2_calculate_means_rbe.py weight_function.csv
  )

foreach(_script ${radioprotection_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS radioprotection DESTINATION bin)

