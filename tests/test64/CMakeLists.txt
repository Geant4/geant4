cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
project( test64 )

find_package( Geant4 REQUIRED )
include( ${Geant4_USE_FILE} )

GEANT4_EXECUTABLE( test64 test64.cc src/*.cc )

#---Test definition-----------------------------------------
GEANT4_ADD_TEST( test64-build BUILD test64 
                              LABELS Nightly Continuous ) 

# Bertini tests
GEANT4_ADD_TEST( test64-bert
                       COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test64 ${CMAKE_CURRENT_SOURCE_DIR}/bert.in
	               DEPENDS test64-build
  	               LABELS Nightly Continuous
	               ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} )

# Fritiof tests
GEANT4_ADD_TEST( test64-ftfp
                       COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test64 ${CMAKE_CURRENT_SOURCE_DIR}/ftfp.in
	               DEPENDS test64-build
  	               LABELS Nightly Continuous
	               ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} )
