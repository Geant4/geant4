# Copied from Geant4 installation, need patch
# - Find StatTest
# This module tries to find the StatTest application
# Once done this will define
#
#	STATTEST_FOUND    - Application found
#	STATTEST_APP      - Application
#   STATTEST_CMD      - Command line to run the application
#   STATTEST_ADD_TEST - Helper function to define a ctest using
#      			          StatTest (Only if we are using for G4 testing)
#
# Variables used by this module, which can change the default behaviour
# and need to be set before calling find_package
#
#	STATTEST_ROOT_DIR	Root directory to StatTest package
#

if(NOT GEANT4_ENABLE_TESTING)
    #No meaning for this if not testing
    SET(STATTEST_FOUND FALSE)
    return()
endif()

#Search application
#Note that the second suggested path is G4 specific...
find_path(STATTEST_APP_DIR
		NAMES StatTestVersion.py
		PATHS ${STATTEST_ROOT_DIR} ${PROJECT_SOURCE_DIR}/verification/StatTest
		NO_DEFAULT_PATH
    )

#If we didn't find it there, fall back to some standard search paths
find_path(STATTEST_APP_DIR NAMES StatTestVersion.py)

####
#### Run-time dependencies...
####
# - Check if ROOT is available
find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
  MESSAGE(STATUS "StatTest: ROOT package not found --> StatTest package disabled")
  set(STATTEST_FOUND FALSE)
  set(_root_isok FALSE)
else()
  set(_root_isok TRUE)
  # - Check if python interpreter is correct version
  # (the one compatible with root)
  find_package(Python ${ROOT_PYTHONVER} QUIET COMPONENTS Interpreter)
  if(NOT Python_Interpreter_FOUND)
     MESSAGE(STATUS "StatTest: Python interpreter of version ${ROOT_PYTHONVER} not found --> StatTest package disabled")
     set(STATTEST_FOUND FALSE)
     set(_python_isok FALSE)
  else()
     set(STATTEST_FOUND TRUE)
     set(_python_isok TRUE)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    StatTest
    DEFAULT_MSG
    STATTEST_APP_DIR
    _root_isok
    _python_isok
  )


if(STATTEST_FOUND)
  set(STATTEST_APP ${STATTEST_APP_DIR}/runtests.py)
  set(STATTEST_CMD ${PYTHON_EXECUTABLE} ${STATTEST_APP})

  # Let's create a function that helps in building tests for
  # regression testing
  # function STATTEST_ADD_TEST(<name>
  #	   		                     G4TEST testname
  #                            CONFIG conffile
  #                            INPUT inputfile
  #			         [DEBUG]
  # 				 [TEXT]  #Input is text file instead of root
  #                            [REFERENCE reference]
  #                            [LABELS label1 label2 ...]
  #                            [IMG filename])
  function(STATTEST_ADD_TEST stattest)
        CMAKE_PARSE_ARGUMENTS(ARG "DEBUG;TEXT"
        "CONFIG;INPUT;REFERENCE;G4TEST;IMG"
        "LABELS" ${ARGN}
    )
    #check mandatory arguments
    list(LENGTH ARG_CONFIG _len)
    if(_len LESS 1)
        message(FATAL_ERROR "STATTEST_ADD_TEST: conffile is mandatory")
    endif()

    list(LENGTH ARG_INPUT _len)
    if(_len LESS 1)
        message(FATAL_ERROR "STATTEST_ADD_TEST: inputfile is mandatory")
    endif()

    list(LENGTH ARG_G4TEST _len)
    if(_len LESS 1)
        message(FATAL_ERROR "STATTEST_ADD_TEST: testname is mandatory")
    endif()

    #Set basic command line
    if(ARG_IMG)
        set(_command ${STATTEST_CMD} -g ${ARG_IMG})
    else()
        set(_command ${STATTEST_CMD})
    endif()

    if(ARG_TEXT)
	    set(_command ${_command} "-T")
    endif()

    #Mandatory parameters
    set(_command ${_command} ${ARG_CONFIG} ${ARG_INPUT})
    if(ARG_REFERENCE)
        set(_command ${_command} ${ARG_REFERENCE})
    endif()

    if(ARG_LABELS)
        set(_labels ${ARG_LABELS})
    else()
        set(_labels "")
    endif()

    include(G4CTest)

     #Now build G4 test
    GEANT4_ADD_TEST(${stattest}
        COMMAND ${_command}
        DEPENDS ${ARG_G4TEST}
        LABELS ${_labels}
        )
  endfunction()
endif(STATTEST_FOUND)

mark_as_advanced(STATTEST_APP_DIR)
