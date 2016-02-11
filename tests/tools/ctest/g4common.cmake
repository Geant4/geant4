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

#---Set the CTEST SITE according to the environment-------------------------
if("$ENV{CTEST_SITE}" STREQUAL "")
  set(CTEST_SITE "${host}")
else()
  set(CTEST_SITE "$ENV{CTEST_SITE}")
  message( "Running build and test on ${host}" )
endif()

#---------------------------------------------------------------------------

if(WIN32)
  if(tag MATCHES x64)
    set(win64 " Win64")
  endif()
  # be4 adding a new generator, make sure that cmake knows about this....
  if(tag MATCHES vc14)
    set(CTEST_CMAKE_GENERATOR "Visual Studio 14 2015${win64}")
  elseif(tag MATCHES vc12)
    set(CTEST_CMAKE_GENERATOR "Visual Studio 12${win64}")
  elseif(tag MATCHES vc11)
    set(CTEST_CMAKE_GENERATOR "Visual Studio 11${win64}")
  elseif(tag MATCHES vc10)
    set(CTEST_CMAKE_GENERATOR "Visual Studio 10${win64}")
  elseif(tag MATCHES vc9)
    set(CTEST_CMAKE_GENERATOR "Visual Studio 9 2008${win64}")
  else()
    set(CTEST_CMAKE_GENERATOR "NMake Makefiles")
  endif()
else()
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
  set(CTEST_BUILD_COMMAND "make -s -i -j${ncpu}") 
endif()

list(APPEND CMAKE_CONFIGURATION_TYPES Release Debug RelWithDebInfo MinSizeRel TestRelease Maintainer)

if("$ENV{BUILDTYPE}" STREQUAL "" OR "$ENV{BUILDTYPE}" STREQUAL "RelWithDebInfo")
  set(CTEST_BUILD_CONFIGURATION "RelWithDebInfo")
  set(CTEST_BUILD_NAME ${tag})
elseif("${CMAKE_CONFIGURATION_TYPES}" MATCHES "$ENV{BUILDTYPE}")
  set(CTEST_BUILD_CONFIGURATION "$ENV{BUILDTYPE}")
  set(CTEST_BUILD_NAME ${tag}-$ENV{BUILDTYPE})
else()
  set(CTEST_BUILD_CONFIGURATION "RelWithDebInfo")
  set(CTEST_BUILD_NAME ${tag}-$ENV{BUILDTYPE})
endif()
if(NOT $ENV{VERSION} STREQUAL "g4tags-dev")
  set(CTEST_BUILD_NAME $ENV{VERSION}-${CTEST_BUILD_NAME})
endif()

if (NOT "$ENV{THREAD}" STREQUAL "")
  set (CTEST_BUILD_NAME ${CTEST_BUILD_NAME}-$ENV{THREAD})
endif()

if (NOT "$ENV{BUILDOPTIONS}" STREQUAL "")
  set (CTEST_BUILD_NAME ${CTEST_BUILD_NAME}-$ENV{BUILDOPTIONS})
endif()

set(CTEST_CONFIGURATION_TYPE "${CTEST_BUILD_CONFIGURATION}")

#---CDash settings----------------------------------------------------------
set(CTEST_PROJECT_NAME "Geant4")
set(CTEST_NIGHTLY_START_TIME "00:00:00 CET")
set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "cdash.cern.ch")
set(CTEST_DROP_LOCATION "/submit.php?project=Geant4")
set(CTEST_DROP_SITE_CDASH TRUE)

#---Custom CTest settings---------------------------------------------------
set(CTEST_CUSTOM_TESTS_IGNORE test19 test29 test39 test49 test47
                              example-ext-biasing-b02 example-adv-eRosita)
set(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE "100000")
set(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE "10000")
if(WIN32)
  set(CTEST_CUSTOM_TESTS_IGNORE ${CTEST_CUSTOM_TESTS_IGNORE} example-ext-geometry-olap)
  set(CTEST_CUSTOM_WARNING_EXCEPTION ${CTEST_CUSTOM_WARNING_EXCEPTION}
        "Ranlux64Engine.+: warning C4293:"
	      "SystemOfUnits.+: warning C4005: 'pascal'"
        ": warning LNK4221:")
else()
  set(CTEST_CUSTOM_WARNING_EXCEPTION ${CTEST_CUSTOM_WARNING_EXCEPTION}
        "warning: ignoring return value of"
        "clang: warning: argument unused" 
        "include/xercesc/util/regx/Token.hpp"
		  "note: variable tracking size limit exceeded with -fvar-tracking-assignments")
endif()
#	     "warning: declaration of 'tokType' shadows a member of"

#message( "CTEST_TIMEOUT =x$ENV{CTEST_TIMEOUT}x")
if("$ENV{CTEST_TIMEOUT}" STREQUAL "")
  # use default value
   set(CTEST_TEST_TIMEOUT 1500)
else()
   message( "Setting timelimit from environment to $ENV{CTEST_TIMEOUT}" )
   set(CTEST_TEST_TIMEOUT $ENV{CTEST_TIMEOUT})
endif()

if("$ENV{VERSION}" STREQUAL "g4tags-dev")
  set(CTEST_NOTES_FILES ${CTEST_SOURCE_DIRECTORY}/gettags.txt)
endif()
set(CTEST_NOTES_FILES ${CTEST_NOTES_FILES} ${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME})

if(EXISTS ${pwd}/jenkins.err)
  set(CTEST_NOTES_FILES ${CTEST_NOTES_FILES} ${pwd}/jenkins.err)
endif()
if(EXISTS ${pwd}/jenkins.log)
  set(CTEST_NOTES_FILES ${CTEST_NOTES_FILES} ${pwd}/jenkins.log)
endif()
  

set(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY_ONCE 1)
set(CTEST_CONFIG_OPTIONS -DGEANT4_ENABLE_TESTING=ON
                         -DGEANT4_BUILD_EXAMPLES=OFF
                         -DGEANT4_INSTALL_DATA=ON
                         -DGEANT4_USE_GDML=ON
                         -DXERCESC_ROOT_DIR=$ENV{XERCESC_ROOT_DIR}
                         $ENV{G4_XOPTS})
if(WIN32)
  set(CTEST_UPDATE_COMMAND ${CTEST_SCRIPT_DIRECTORY}/g4tagsvn.bat)
else()
  set(CTEST_UPDATE_COMMAND ${CTEST_SCRIPT_DIRECTORY}/g4tagsvn.py)
endif()
set(CTEST_UPDATE_TYPE SVN)
set(CTEST_UPDATE_OPTIONS "-c $ENV{VERSION} -d ${CTEST_SOURCE_DIRECTORY}")

#---Set Runtime environment-------------------------------------------------
if(WIN32) 
  if(NOT CTEST_CMAKE_GENERATOR MATCHES Makefiles)
    set(_cfg /${CTEST_BUILD_CONFIGURATION})
  endif()
#  set(ENV{PATH} "${CTEST_BINARY_DIRECTORY}/outputs/runtime${_cfg};$ENV{PATH}")
  set(ENV{PATH} "${CTEST_BINARY_DIRECTORY}/BuildProducts/${_cfg}/bin;$ENV{PATH}")
endif()
