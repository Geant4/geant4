#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test14)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
GEANT4_EXECUTABLE(test14 test14.cc src/*.cc)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build test14. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(INPUTS 
     test14_lowE.large_N.in
	  test14_penelope.large_N.in
	  test14_polarised.large_N.in
	  )
	  
foreach(_input ${INPUTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_input}
    ${PROJECT_BINARY_DIR}/${_input}
    COPYONLY
    )
endforeach()


# add the tests to build and run the executable

GEANT4_ADD_TEST(test14-build BUILD test14)

GEANT4_ADD_TEST(test14-lowE 
                  COMMAND test14 test14_lowE.large_N.in
						ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
						DEPENDS test14-build
					)
					
GEANT4_ADD_TEST(test14-penelope
                  COMMAND test14 test14_penelope.large_N.in
						ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
						DEPENDS test14-build
					)
					
GEANT4_ADD_TEST(test14-polarised
                  COMMAND test14 test14_polarised.large_N.in
						ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
						DEPENDS test14-build
					)
