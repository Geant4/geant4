#
# EXAMPLE USAGE WITH G4.9.6.ref02:
#
# mkdir /tests/test75-build
# cd /tests/test75-build
# cmake -DGeant4_SOURCE_DIR=$G4INSTALL \
#       -DGeant4_DIR=$G4LIB/Geant4-9.6.0 \
#       ../test75
# make
# make install
#
# NOTE-1: it explicitly requires out-of-source build
# NOTE-2: everything (incl. configs, Root macros, etc.) will be copied to test75-build/bin
#
# Note-3: somehow, passing in ROOT location is NOT necessary
# Note-4: XERCESCROOT (if needed) gets setup by the like of "setup xerces_c v3_1_1..."
#
cmake_minimum_required(VERSION 2.6.4 FATAL_ERROR)

#----Included in Geant4------------------------------------------------------------------
if(CMAKE_PROJECT_NAME STREQUAL Geant4)
  find_package(Geant4 REQUIRED)
  include(${Geant4_USE_FILE})

  find_package(ROOT QUIET)
  if(NOT ROOT_FOUND)
    message(STATUS "G4 TESTS: ROOT package not found. --> test75 disabled")
    return()
  endif()

  GEANT4_ADD_TEST(test75-build BUILD test75)
  GEANT4_ADD_TEST(test75-gam300 COMMAND test75 ${CMAKE_CURRENT_SOURCE_DIR}/ctest75-gam300-Cu.config
                         DEPENDS test75-build
                         ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

else()
#----Standalone--------------------------------------------------------------------------
  project(test75)

  if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "test75 requires an out-of-source build.")
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
    message(STATUS "G4 TESTS: ROOT package not found. --> test75 disabled")
    return()
  endif()
    
  #----------------------------------------------------------------------------
  # Locate sources and headers for this project
  # NB: headers are included so they will show up in IDEs
  #
  file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
  file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
  
  #----------------------------------------------------------------------------
  # Add the executable, and link it to the Geant4 libraries
  #
  add_executable(test75 test75.cc ${sources} ${headers})
  
  include_directories(${PROJECT_SOURCE_DIR}/include ${INCLUDE_DIRECTORIES} ${EXTRA_INCLUDE_DIRECTORIES} ${ROOT_INCLUDE_DIR})

  target_link_libraries(test75 ${EXTRA_LIBRARIES} ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${XercesC_LIBRARIES})
    
  
  #----------------------------------------------------------------------------
  # Copy all scripts, input files to the build directory, i.e. the directory in which we
  # build test75.
  #
  #  Note that it is not recommend to use GLOB to collect a list of
  #  source files etc... from your source tree. If no CMakeLists.txt file
  #  changes when a source is added or removed then the generated build
  #  system cannot know when to ask CMake to regenerate.
  
  file(GLOB TEST75_FILES RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/test75.gamma
                                                        ${PROJECT_SOURCE_DIR}/*.C
                                                        ${PROJECT_SOURCE_DIR}/GammaData/*.C
                                                        README
                                                        )
  
  foreach(_file ${TEST75_FILES})
    configure_file(
      ${PROJECT_SOURCE_DIR}/${_file} ${PROJECT_BINARY_DIR}/${_file} COPYONLY
      )
  endforeach()
  
  #----------------------------------------------------------------------------
  # Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
  #
  set(CMAKE_INSTALL_PREFIX ./)
  install(TARGETS test75 DESTINATION bin)
  
endif()
