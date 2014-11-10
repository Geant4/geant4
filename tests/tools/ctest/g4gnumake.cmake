cmake_minimum_required(VERSION 2.6)

include(${CTEST_SCRIPT_DIRECTORY}/g4macros.cmake)

GET_PWD(pwd)
GET_HOST(host)
GET_NCPUS(ncpu)
GET_CONFIGURATION_TAG(tag)

set(MODEL $ENV{MODE})
#---Set the CTEST SITE according to the environment-------------------------
if("$ENV{CTEST_SITE}" STREQUAL "")
  set(CTEST_SITE "${host}")
else()
  set(CTEST_SITE "$ENV{CTEST_SITE}")
  message( "Running build and test on ${host}" )
endif()

set(CTEST_BUILD_NAME ${tag}-gnumake)
if(NOT $ENV{VERSION} STREQUAL "g4tags-dev")
  set(CTEST_BUILD_NAME $ENV{VERSION}-${CTEST_BUILD_NAME})
endif()


if(DEFINED ENV{WITH_MEMCHECK})
  set(WITH_MEMCHECK TRUE)
  find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
  set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--track-origins=yes")
endif()

set(CTEST_SOURCE_DIRECTORY "$ENV{SOURCE}")
set(CTEST_BINARY_DIRECTORY "$ENV{BINARY}")
set(CTEST_DASHBOARD_ROOT "$ENV{BINARY}/dashboard")

#--Customize to GNU make build commands--------------------------------------
#set(CTEST_CHECKOUT_COMMAND ${CTEST_SCRIPT_DIRECTORY}/g4tagsvn.py checkout -c $ENV{VERSION} -d $ENV{SOURCE} -q)
set(CTEST_UPDATE_COMMAND ${CTEST_SCRIPT_DIRECTORY}/g4tagsvn.py)
set(CTEST_UPDATE_TYPE SVN)
set(CTEST_UPDATE_OPTIONS "-c $ENV{VERSION} -d $ENV{SOURCE}")

set(CTEST_CONFIGURE_COMMAND "echo configuration not needed")
set(CTEST_BUILD_COMMAND "make -s -i -j${ncpu}") 
 
#---CDash settings----------------------------------------------------------
set(CTEST_PROJECT_NAME "Geant4")
set(CTEST_NIGHTLY_START_TIME "00:00:00 CET")
set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "cdash.cern.ch")
set(CTEST_DROP_LOCATION "/submit.php?project=Geant4")
set(CTEST_DROP_SITE_CDASH TRUE)

#---Addional CTest settings-------------------------------------------------
#set(CTEST_UPDATE_OPTIONS "${CTEST_UPDATE_OPTIONS} -p")  # Add proposed tags
set(CTEST_TEST_TIMEOUT 1500)
set(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE "100000")
set(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE "10000")

set($ENV{LC_MESSAGES} "en_EN")
set(CTEST_CUSTOM_WARNING_EXCEPTION ${CTEST_CUSTOM_WARNING_EXCEPTION}
		  "note: variable tracking size limit exceeded with -fvar-tracking-assignments")

#---Configure tests. Some of them require some files to be copied-----------
configure_file(${CTEST_SOURCE_DIRECTORY}/tests/ctests/GNUMakeTestfile.cmake ${CTEST_BINARY_DIRECTORY}/CTestTestfile.cmake COPYONLY)
configure_file(${CTEST_SOURCE_DIRECTORY}/examples/extended/runAndEvent/RE05/pythia_event.data ${CTEST_BINARY_DIRECTORY}/pythia_event.data)

#---CTest commands----------------------------------------------------------
ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
ctest_start(${MODEL} TRACK ${MODEL})
ctest_update(SOURCE ${CTEST_SOURCE_DIRECTORY} )
ctest_configure(BUILD ${CTEST_BINARY_DIRECTORY} )
ctest_build(BUILD ${CTEST_SOURCE_DIRECTORY}/source)
ctest_test(PARALLEL_LEVEL ${ncpu})

if (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
   set(CTEST_TEST_TIMEOUT 3500)
   ctest_memcheck(PARALLEL_LEVEL ${ncpu} 
			  INCLUDE "example-(bas|ext)"
			  EXCLUDE "-build$")
endif (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
ctest_submit()
