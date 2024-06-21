#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(HepMCEx02)

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
# Find Pythia6 (optional package)
#
find_package(Pythia6 QUIET)
if(Pythia6_FOUND)
  message(STATUS "G4 Examples: Pythia6 found. --> HepMCEx02 example with Pythia6 enabled.") 
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${HEPMC_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(HepMCEx02 HepMCEx02.cc ${sources} ${headers})
target_compile_definitions(example_MyPythia PRIVATE $<$<BOOL:${Pythia6_FOUND}>:G4LIB_USE_PYTHIA>)
#target_link_libraries(HepMCEx02 ${Geant4_LIBRARIES}
#                      ${HEPMC_LIBRARIES} ${HEPMC_FIO_LIBRARIES}
#                      $<$<BOOL:${Pythia6_FOUND}>:Pythia6::Pythia6>
#                      gfortran)
target_link_libraries(HepMCEx02 ${Geant4_LIBRARIES}
                      ${HEPMC_LIBRARIES} ${HEPMC_FIO_LIBRARIES}
                      $<$<BOOL:${Pythia6_FOUND}>:Pythia6::Pythia6>)

# if pythia is compiled with g77, link with -lg2c instead.
#target_link_libraries(HepMCEx02 ${Geant4_LIBRARIES}
#                      ${HEPMC_LIBRARIES} ${HEPMC_FIO_LIBRARIES}
#                      $<$<BOOL:${Pythia6_FOUND}>:Pythia6::Pythia6>
#                      g2c)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build HepMCEx02. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(HepMCEx02_SCRIPTS
    hepmc_ascii.in hepmc_ascii.out hepmc_pygen.in hepmc_pygen.out pyh4mu.mac pyset.mac vis.mac
    data/example_MyPythia.dat
  )

foreach(_script ${HepMCEx02_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS HepMCEx02 DESTINATION bin)

