cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test30)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test30 test30.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test30 COMMAND test30 ${CMAKE_CURRENT_SOURCE_DIR}/test30.in
                       BUILD test30
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

#---Physics Validaiton-------------------------------------------------------------------------------
# Tests executed only in ctest/cdash PhysicsChecks group
# Testing for PhysicsChecks is done in three steps:
# 1- A test trying to build the application
# 2- A simulation job
# 3- A step in which the output of the simulation is checked against a reference

#Step1 Build executable test
GEANT4_ADD_TEST(test30-build BUILD test30 LABELS PhysicsChecks)

#Step 2 Perform simulation test: runs all models together
GEANT4_ADD_TEST(test30-Bi 
			  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test30 ${CMAKE_CURRENT_SOURCE_DIR}/physval/macro.in 
			  DEPENDS test30-build
			  LABELS PhysicsChecks
			  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

#Step 3 Execute regression testing only if StatTest program is available
#Each "model" is a separate test
find_package(StatTest QUIET)
if(STATTEST_FOUND)
	foreach( _model bertini binary bertini_preco)
		 STATTEST_ADD_TEST( test30-Bi-${_model}-checkOutput
				G4TEST test30-Bi
				CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/physval/testconf.qa
				INPUT ${_model}.root
				REFERENCE ${CMAKE_CURRENT_SOURCE_DIR}/physval/${_model}-ref.root
				IMG ${_model}.pdf
				LABELS PhysicsChecks)
	endforeach()
endif()
