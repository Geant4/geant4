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
# Find ROOT variables if the variable G4ANALYSIS_USE_ROOT is set
#
find_package(ROOT QUIET)


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
if(ROOT_FOUND)
add_definitions(-DG4ANALYSIS_USE_ROOT)
include(${Geant4_USE_FILE})
include_directories(${ROOT_INCLUDE_DIR} ${Geant4_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/include)
else()
include(${Geant4_USE_FILE})
<<<<<<< HEAD
include_directories(${PROJECT_SOURCE_DIR}/include)
endif()

=======

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
# find_package(ROOT REQUIRED)
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

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

if(ROOT_FOUND)
target_link_libraries(hadrontherapy ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
else()
target_link_libraries(hadrontherapy ${Geant4_LIBRARIES})
<<<<<<< HEAD
endif()
=======
# ${ROOT_LIBRARIES} )
file(GLOB_RECURSE dataFiles RELATIVE ${PROJECT_SOURCE_DIR} data/rbe/*.*)
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build hadrontherapy. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
<<<<<<< HEAD
set(EXAMPLEhadrontherapy_SCRIPTS
  clean.sh
 batch.mac
 macro/GUIPersonalisation.mac
  macro/carbon_beamline.mac
  macro/detectorGeometry.mac
  macro/hadron_therapy.mac
  macro/manySlabsPhantom.mac
  macro/modulator.loop
  macro/modulatorMacro.mac
  macro/oneSlabPhantom.mac
  macro/stoppingPowers.mac
  macro/visualisationMacro.mac
  macro/voxellizedPhantom.mac
  macro/laserDrivenBeamline.mac
  macro/defaultMacroWithReferencePhysicsList.mac
  macro/defaultMacro.mac
  macro/Taranis.mac
  macro/inputFC.mac

 field/ESSMagneticField.TABLE

# Faraday Cup
 field/ElectricFieldFC-600V.TABLE

#QUADRUPOLE 40mm
 field/Quad40MagneticField.TABLE
=======
file(GLOB macros RELATIVE ${PROJECT_SOURCE_DIR} macro/*.mac  macro/*.loop)
file(GLOB modulators RELATIVE ${PROJECT_SOURCE_DIR} Modulators/*.txt)
file(GLOB fields RELATIVE ${PROJECT_SOURCE_DIR} field/*.TABLE)
file(GLOB_RECURSE data RELATIVE ${PROJECT_SOURCE_DIR} data/*)
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

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

