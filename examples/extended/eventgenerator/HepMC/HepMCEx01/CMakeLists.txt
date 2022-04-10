#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(HepMCEx01 LANGUAGES C CXX Fortran )

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
# Find HepMC (required package)
#
find_package(HepMC REQUIRED)
#----------------------------------------------------------------------------
# Find LHAPDF (required package)
#
find_package(LHAPDF REQUIRED)
#----------------------------------------------------------------------------
# Find Pythia6 (optional package)
#
find_package(Pythia6 QUIET)
if(PYTHIA6_FOUND)
  message(STATUS "G4 Examples: Pythia6 found. --> HepMCEx01 example with Pythia6 enabled.") 
  add_definitions(-DG4LIB_USE_PYTHIA)
    enable_language(Fortran)
  set (generator_sources 
   ${PROJECT_SOURCE_DIR}/external/upevnt.F
   ${PROJECT_SOURCE_DIR}/external/upinit.F
   ${PROJECT_SOURCE_DIR}/external/upveto.F
   )
else()
  set(PYTHIA6_LIBRARIES "")  
  set (generator_sources "")
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(HepMCEx01 HepMCEx01.cc ${sources} ${headers} ${generator_sources})
target_link_libraries(HepMCEx01 ${Geant4_LIBRARIES}
                      ${HEPMC_LIBRARIES} ${HEPMC_FIO_LIBRARIES}
                      ${PYTHIA6_LIBRARIES} ${LHAPDF_LIBRARIES})
target_include_directories(HepMCEx01 PRIVATE ${PROJECT_SOURCE_DIR}/include 
                                             ${Geant4_INCLUDE_DIR}
                                             ${HEPMC_INCLUDE_DIR})
if(PYTHIA6_FOUND)
set_target_properties(HepMCEx01 PROPERTIES LINKER_LANGUAGE Fortran)
endif()
#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build HepMCEx01. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(HepMCEx01_SCRIPTS
    hepmc_ascii.in hepmc_ascii.out hepmc_pygen.in hepmc_pygen.out pyh4mu.mac pyset.mac vis.mac
    data/example_MyPythia.dat
  )

foreach(_script ${HepMCEx01_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS HepMCEx01 DESTINATION bin)

