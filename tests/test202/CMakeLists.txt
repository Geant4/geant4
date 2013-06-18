cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test202)
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test202 test202.cc src/*.cc)

configure_file(testVisSceneAdd.mac testVisSceneAdd.mac COPYONLY) # test expects file to be local
configure_file(testVisSet.mac testVisSet.mac COPYONLY) # test expects file to be local
configure_file(vis.mac vis.mac COPYONLY) # test expects file to be local
configure_file(run0.mac run0.mac COPYONLY) # this file local too

#---Test definitions----------------------------------------------------------------------------------

#build example:
GEANT4_ADD_TEST(test202-build 
     BUILD test202
     SOURCE_DIR ${CMAKE_SOURCE_DIR}/tests/test202
     PROJECT test202)
