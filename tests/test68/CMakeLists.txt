cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
project( test68 )

find_package( Geant4 REQUIRED )
include( ${Geant4_USE_FILE} )

find_package( PythonInterp QUIET )
if( NOT PYTHONINTERP_FOUND )
  message( STATUS "Python interpreter not found !" )
  return()
endif()

GEANT4_EXECUTABLE( test68 test68.cc src/*.cc )

#---Test definition-----------------------------------------
GEANT4_ADD_TEST( test68-build BUILD test68 
                              LABELS Experimental ) 

#message( ${CMAKE_CURRENT_SOURCE_DIR} " " ${CMAKE_CURRENT_BINARY_DIR} )

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/start.rndm ${CMAKE_CURRENT_BINARY_DIR}/start.rndm COPYONLY )

GEANT4_ADD_TEST( test68-reproducibility 
		 COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/reproducibility.py 
		 DEPENDS test68-build
		 LABELS Experimental
		 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} )
