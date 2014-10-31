#
# This is a crude yet working version; tested w/g4.9.5.p01
# Initially modeled after examples/extended/persistency/P01/CMakeLists.txt
# Modified to include mucapture (default is OFF) and enable/force an out-of-source build
#
# USAGE: cmake -DGeant4_DIR=location-of-geant4-binaries \
#              -DGeant4_SOURCE_DIR=location-of-geant4-source\
#              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
#              -DUSE_MUCAPTURE=ON \
#              -DMuCapture_SOURCE_DIR=location-of-mucapture-source \
#              -DXERCESC_ROOT_DIR=$XERCESCROOT \
#              location_of_test48_source_itself
# 
#
# EXAMPLE USAGE WITH G4.9.6.b01:
#
# mkdir /tests/test48-build
# cd /tests/test48-build
# cmake -DGeant4_SOURCE_DIR=$G4INSTALL \
#       -DGeant4_DIR=$G4LIB/Geant4-9.6.0 \
#       -DXERCESC_ROOT_DIR=$XERCESCROOT \
#       ../test48
# make
# make install
#
# NOTE-1: it explicitly requires out-of-source build
# NOTE-2: everything (incl. configs, Root macros, etc.) will be copied to test48-build/bin
#
# Note-3: somehow, passing in ROOT location is NOT necessary
# Note-4: XERCESCROOT gets setup by the like of "setup xerces_c v3_1_1..."
#
cmake_minimum_required(VERSION 2.6.4 FATAL_ERROR)

#----Included in Geant4------------------------------------------------------------------
if(CMAKE_PROJECT_NAME STREQUAL Geant4)
  find_package(Geant4 REQUIRED)
  include(${Geant4_USE_FILE})

  find_package(ROOT QUIET)
  if(NOT ROOT_FOUND)
    message(STATUS "G4 TESTS: ROOT package not found. --> test48 disabled")
    return()
  endif()

  GEANT4_ADD_TEST(test48-build BUILD test48)
  GEANT4_ADD_TEST(test48-piminus COMMAND test48 ${CMAKE_CURRENT_SOURCE_DIR}/ctest48-picapture-Pb.config
                         DEPENDS test48-build
                         ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(test48-muminus COMMAND test48 ${CMAKE_CURRENT_SOURCE_DIR}/ctest48-mucapture-Si.config
                         DEPENDS test48-build
                         ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

else()
#----Standalone--------------------------------------------------------------------------
  project(test48)

  if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "test48 requires an out-of-source build.")
    message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
    message(STATUS "CMakeCache.txt")
    message(STATUS "CMakeFiles")
    message(STATUS "Once these files are removed, create a separate directory")
    message(STATUS "and run CMake from there")
    message(FATAL_ERROR "in-source build detected")
  endif()

  find_package(Geant4 REQUIRED)
  include(${Geant4_USE_FILE})

  set(Geant4_SOURCE_DIR  "$ENV{G4INSTALL}" CACHE PATH "Default location of Geant4 source")
  if (Geant4_SOURCE_DIR)
    set(CMAKE_MODULE_PATH ${Geant4_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})
  endif (Geant4_SOURCE_DIR)

  find_package(ROOT QUIET)
  if(NOT ROOT_FOUND)
    message(STATUS "G4 TESTS: ROOT package not found. --> test48 disabled")
    return()
  endif()
  #find_package(XercesC REQUIRED)

  #set(CMAKE_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/include")  
  #message("CMAKE_INCLUDE_PATH: ${CMAKE_INCLUDE_PATH}")
  #INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/include)

  option (USE_MUCAPTURE "Use MuCapture" OFF)  # turn it off by deafult
  
  # we tell the system where the mucapture source is
  set(MuCapture_SOURCE_DIR ../mucapture CACHE PATH "Default location of mucapture source")
  
  if (USE_MUCAPTURE)
    # we will build mucapture in a subdirectory of the main build 
    add_subdirectory ( ${MuCapture_SOURCE_DIR} ${CMAKE_BINARY_DIR}/mucapture EXCLUDE_FROM_ALL)
    # we assume MuCapture_INCLUDE_DIRECTORIES is defined by the the script in MuCapture_SOURCE_DIR
    #  message("MuCapture_INCLUDE_DIRECTORIES:  ${MuCapture_INCLUDE_DIRECTORIES}")
  
    if(NOT DEFINED MuCapture_INCLUDE_DIRECTORIES)
      message(FATAL_ERROR "MuCapture_INCLUDE_DIRECTORIES not defined")
    endif(NOT DEFINED MuCapture_INCLUDE_DIRECTORIES)
  
    set (EXTRA_INCLUDE_DIRECTORIES "${MuCapture_INCLUDE_DIRECTORIES}")
    set (EXTRA_LIBRARIES mucapture)
  
    add_definitions(-DUSE_MUCAPTURE)
  endif (USE_MUCAPTURE)
    
  #----------------------------------------------------------------------------
  # Locate sources and headers for this project
  # NB: headers are included so they will show up in IDEs
  #
  file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
  file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
  
  #----------------------------------------------------------------------------
  # Add the executable, and link it to the Geant4 libraries
  #
  add_executable(test48 test48.cc ${sources} ${headers})
  
  include_directories(${PROJECT_SOURCE_DIR}/include ${INCLUDE_DIRECTORIES} ${EXTRA_INCLUDE_DIRECTORIES} ${ROOT_INCLUDE_DIR})

  target_link_libraries(test48 ${EXTRA_LIBRARIES} ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${XercesC_LIBRARIES})
  
  #message("CMAKE_INCLUDE_PATH: ${CMAKE_INCLUDE_PATH}")
  #message("INCLUDE_DIRECTORIES: ${INCLUDE_DIRECTORIES}")
  
  
  #----------------------------------------------------------------------------
  # Copy all scripts, input files to the build directory, i.e. the directory in which we
  # build test48.
  #
  #  Note that it is not recommend to use GLOB to collect a list of
  #  source files etc... from your source tree. If no CMakeLists.txt file
  #  changes when a source is added or removed then the generated build
  #  system cannot know when to ask CMake to regenerate.
  
  file(GLOB TEST48_FILES RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/test48.*minus* 
                                                        ${PROJECT_SOURCE_DIR}/ctest48-picapture-Pb.config
                                                        ${PROJECT_SOURCE_DIR}/*.C
                                                        ${PROJECT_SOURCE_DIR}/antiproton/*.dat
                                                        ${PROJECT_SOURCE_DIR}/kminus/*.dat
                                                        ${PROJECT_SOURCE_DIR}/muminus/*.dat
                                                        ${PROJECT_SOURCE_DIR}/piminus/*.dat
                                                        test48.antiproton
							${PROJECT_SOURCE_DIR}/batch/*.sh
							${PROJECT_SOURCE_DIR}/g4val-upload-scripts/generate*.sh
                                                        README
                                                        )
  
  foreach(_file ${TEST48_FILES})
    configure_file(
      ${PROJECT_SOURCE_DIR}/${_file} ${PROJECT_BINARY_DIR}/${_file} COPYONLY
      )
  endforeach()
  
  #----------------------------------------------------------------------------
  # Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
  #
  set(CMAKE_INSTALL_PREFIX ./)
  # message("CMAKE_INSTALL_PREFIX now set to: ${CMAKE_INSTALL_PREFIX}")
  install(TARGETS test48 DESTINATION bin)
  
endif()
