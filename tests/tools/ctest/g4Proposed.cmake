cmake_minimum_required(VERSION 2.6)

#---Utility Macros----------------------------------------------------------
include(${CTEST_SCRIPT_DIRECTORY}/g4macros.cmake)

#---Make sure that VERBOSE is OFF to avoid screwing up the build performance
unset(ENV{VERBOSE})

#---General Configuration---------------------------------------------------
GET_PWD(pwd)
GET_HOST(host)
GET_NCPUS(ncpu)
GET_CONFIGURATION_TAG(tag)

#---Set the source and build directory--------------------------------------
if("$ENV{SOURCE}" STREQUAL "")
  set(CTEST_SOURCE_DIRECTORY "${pwd}/source")
else()
  set(CTEST_SOURCE_DIRECTORY "$ENV{SOURCE}")
endif()
if("$ENV{BINARY}" STREQUAL "")
  set(CTEST_BINARY_DIRECTORY "${pwd}/build")
else()
  set(CTEST_BINARY_DIRECTORY "$ENV{BINARY}")
endif()

#---------------------------------------------------------------------------

set(CTEST_SITE "${host}")
if(WIN32)
  if(tag MATCHES vc10)
    set(CTEST_CMAKE_GENERATOR "Visual Studio 10")
  elseif(tag MATCHES vc9)
    set(CTEST_CMAKE_GENERATOR "Visual Studio 9 2008")
  else()
    set(CTEST_CMAKE_GENERATOR "NMake Makefiles")
  endif()
else()
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
  set(CTEST_BUILD_COMMAND "make -s -i -j${ncpu}") 
endif()
set(CTEST_BUILD_CONFIGURATION "RelWithDebInfo")
set(CTEST_CONFIGURATION_TYPE "${CTEST_BUILD_CONFIGURATION}")
set(CTEST_BUILD_NAME ${tag}-proposed)

#---CDash settings----------------------------------------------------------
set(CTEST_PROJECT_NAME "Geant4")
set(CTEST_NIGHTLY_START_TIME "00:00:00 CET")
set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "cdash.cern.ch")
set(CTEST_DROP_LOCATION "/submit.php?project=Geant4")
set(CTEST_DROP_SITE_CDASH TRUE)

#---Custom CTest settings---------------------------------------------------
set(CTEST_CUSTOM_TESTS_IGNORE test19 test29 test39 test49 
                              example-ext-biasing-b02 example-adv-eRosita)
set(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE "20000")
set(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE "10000")
if(WIN32)
  set(CTEST_CUSTOM_TESTS_IGNORE ${CTEST_CUSTOM_TESTS_IGNORE} example-ext-geometry-olap)
endif()

set(CTEST_TEST_TIMEOUT 1200)
set(CTEST_NOTES_FILES  ${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME} ${CTEST_SOURCE_DIRECTORY}/gettags.txt)
set(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY_ONCE 1)
set(CTEST_CONFIG_OPTIONS -DGEANT4_ENABLE_TESTING=ON
                         -DGEANT4_BUILD_EXAMPLES=ON
                         -DGEANT4_INSTALL_DATA=ON
                         -DGEANT4_USE_GDML=ON
                         -DXERCESC_ROOT_DIR=$ENV{XERCESC_ROOT_DIR}
                         $ENV{G4_XOPTS})
set(CTEST_UPDATE_COMMAND ${CTEST_SCRIPT_DIRECTORY}/g4tagsvn.py)
set(CTEST_UPDATE_TYPE SVN)
set(CTEST_UPDATE_OPTIONS "-c $ENV{VERSION} -d ${CTEST_SOURCE_DIRECTORY} -p")

#---Set Runtime environment-------------------------------------------------
if(WIN32) 
  if(NOT CTEST_CMAKE_GENERATOR MATCHES Makefiles)
    set(_cfg /${CTEST_BUILD_CONFIGURATION})
  endif()
  #set(ENV{PATH} "${CTEST_BINARY_DIRECTORY}/outputs/runtime${_cfg};$ENV{PATH}")
  set(ENV{PATH} "${CTEST_BINARY_DIRECTORY}/BuildProducts/${_cfg}/bin;$ENV{PATH}")
endif()

#---CTest commands----------------------------------------------------------
ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
while(${CTEST_ELAPSED_TIME} LESS 72000)
  set(START_TIME ${CTEST_ELAPSED_TIME})
  ctest_start (Continuous)
  #ctest_update(RETURN_VALUE updates)
  file(DOWNLOAD http://lcgapp.cern.ch/spi/cgi-bin/g4tags.py?devline=$ENV{VERSION} ${CTEST_BINARY_DIRECTORY}/newtags.txt STATUS status)
  execute_process(COMMAND diff -I "^#" ${CTEST_BINARY_DIRECTORY}/newtags.txt ${CTEST_SOURCE_DIRECTORY}/gettags.txt
                  RESULT_VARIABLE res OUTPUT_VARIABLE out ERROR_VARIABLE err)
  if(res AND NOT err AND status MATCHES 0)
    #---The list of tags is different. This should trigger an update, re-build and re-test
    #execute_process(COMMAND ${CTEST_SCRIPT_DIRECTORY}/g4checkout.py -c $ENV{VERSION} -d ${CTEST_SOURCE_DIRECTORY} --update --proposed --quiet)
    ctest_update(RETURN_VALUE updates)
    message("Updates: ${updates}")
    ctest_configure(BUILD   ${CTEST_BINARY_DIRECTORY} 
                    SOURCE  ${CTEST_SOURCE_DIRECTORY}
                    OPTIONS "${CTEST_CONFIG_OPTIONS}")
    ctest_build(BUILD ${CTEST_BINARY_DIRECTORY})
    ctest_test(PARALLEL_LEVEL ${ncpu}
               EXCLUDE "largeN$"
               INCLUDE "^test")
    ctest_submit()
  endif()
  ctest_sleep(${START_TIME} 600 ${CTEST_ELAPSED_TIME})
endwhile()
