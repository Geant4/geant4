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
add_executable(hadrontherapy hadrontherapy.cc ${sources} ${headers})
target_link_libraries(hadrontherapy ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build hadrontherapy. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEhadrontherapy_SCRIPTS
  clean.sh
 batch.mac
Modulators/Modulator009_02.txt
Modulators/Modulator010.txt
Modulators/Modulator011.txt
Modulators/ModulatorAM16FR.txt
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
  macro/Trento_parameters.mac

 field/ESSMagneticField.TABLE

# Faraday Cup
 field/ElectricFieldFC-600V.TABLE

#QUADRUPOLE 40mm
 field/Quad40MagneticField.TABLE


#QUADRUPOLE 80mm
 field/Quad80MagneticField.TABLE

  )

foreach(_script ${EXAMPLEhadrontherapy_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(hadrontherapy DEPENDS hadrontherapy)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS hadrontherapy DESTINATION bin)
