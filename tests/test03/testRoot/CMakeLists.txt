# $Id: CMakeLists.txt 68058 2013-03-13 14:47:43Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test03-Root)

#----------------------------------------------------------------------------
# Add path to this test CMake files
#
set(CMAKE_MODULE_PATH
    ${PROJECT_SOURCE_DIR}/../cmake
    ${CMAKE_MODULE_PATH})
    
#----------------------------------------------------------------------------
# Add definitions 
#
add_definitions(-DTEST_ANALYSIS_ROOT)

#----------------------------------------------------------------------------
# Build test for Root output
#

include(BuildTest)

include_directories(${PROJECT_SOURCE_DIR}/../source/include)

GEANT4_EXECUTABLE(test03-Root  ../source/test03.cc ../source/src/*.cc)

GEANT4_ADD_TEST(test03-Root-build BUILD test03-Root)

GEANT4_ADD_TEST(test03-Root-write 
                COMMAND test03-Root  -m run2.mac -read off
                ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                DEPENDS test03-Root-build
                )
						
GEANT4_ADD_TEST(test03-Root-read 
                COMMAND test03-Root  -m run2.mac -read on
                ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                DEPENDS test03-Root-write
                )       
