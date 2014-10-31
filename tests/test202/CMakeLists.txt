cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test202)
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test202 test202.cc src/*.cc)

configure_file(testVisSceneAdd.mac testVisSceneAdd.mac COPYONLY) # test expects file to be local
configure_file(testVisSet.mac testVisSet.mac COPYONLY) # test expects file to be local
configure_file(vis.mac vis.mac COPYONLY) # test expects file to be local
configure_file(run0.mac run0.mac COPYONLY) # this file local too
configure_file(OGLIX.mac OGLIX.mac COPYONLY) # this file local too
configure_file(xvfb.sh xvfb.sh COPYONLY) # this file local too

# OGLIX-reference directory
configure_file(OGLIX-reference/visTest-DrawVolume_0.eps OGLIX-reference/visTest-DrawVolume_0.eps  COPYONLY) # this file local too
configure_file(OGLIX-reference/visTest-VisAdd_0.eps OGLIX-reference/visTest-VisAdd_0.eps  COPYONLY) # this file local too
configure_file(OGLIX-reference/visTest-VisSet_0.eps OGLIX-reference/visTest-VisSet_0.eps  COPYONLY) # this file local too
configure_file(OGLIX-reference/visTest-BeamOn_0.eps OGLIX-reference/visTest-BeamOn_0.eps  COPYONLY) # this file local too

#---Test definitions----------------------------------------------------------------------------

#build example:
GEANT4_ADD_TEST(test202-build 
     BUILD test202
     SOURCE_DIR ${CMAKE_SOURCE_DIR}/tests/test202
     PROJECT test202)

find_package(Xvfb QUIET)

if(XVFB_FOUND)
  message(STATUS "G4 TESTS: found Xvfb --> run test202")

  GEANT4_ADD_TEST(test202-OGLIX COMMAND xvfb.sh ${CMAKE_CURRENT_BINARY_DIR}/test202 OGLIX
                              ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                              DEPENDS test202-build
			      TIMEOUT 1500
			      ERROR expected_errors.txt)

    
endif()
