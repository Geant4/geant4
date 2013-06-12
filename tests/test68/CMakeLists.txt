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
GEANT4_ADD_TEST( test68-build BUILD test68 LABELS Nightly Continuous ) 

#message( ${CMAKE_CURRENT_SOURCE_DIR} " " ${CMAKE_CURRENT_BINARY_DIR} )

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/start.rndm ${CMAKE_CURRENT_BINARY_DIR}/start.rndm COPYONLY )

GEANT4_ADD_TEST( test68-reproducibility-k0lFeSci 
		 COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/reproducibility_k0lFeSci.py $<CONFIGURATION> 
                 TIMEOUT 3600
		 DEPENDS test68-build
		 LABELS Nightly Continuous
		 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} )

GEANT4_ADD_TEST( test68-reproducibility-pimCuLAr 
		 COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/reproducibility_pimCuLAr.py $<CONFIGURATION> 
                 TIMEOUT 3600
		 DEPENDS test68-build
		 LABELS Nightly Continuous
		 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} )

GEANT4_ADD_TEST( test68-reproducibility-nPbSci 
		 COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/reproducibility_nPbSci.py $<CONFIGURATION> 
                 TIMEOUT 3600
		 DEPENDS test68-build
		 LABELS Nightly Continuous
		 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} )

GEANT4_ADD_TEST( test68-reproducibility-pWLAr 
		 COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/reproducibility_pWLAr.py $<CONFIGURATION> 
                 TIMEOUT 3600
		 DEPENDS test68-build
		 LABELS Nightly Continuous
		 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} )

GEANT4_ADD_TEST( test68-reproducibility-kmPbWO4 
		 COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/reproducibility_kmPbWO4.py $<CONFIGURATION> 
                 TIMEOUT 3600
		 DEPENDS test68-build
		 LABELS Nightly Continuous
		 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} )
