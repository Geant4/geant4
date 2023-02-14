# - Internal Use file for Geant4
# This file is designed for inclusion by the UseGeant4.cmake file, but
# only in the build tree. It contains functions and macros that are only
# intended to help in building the Geant4 tests.
#
# IT SHOULD NOT BE INSTALLED!

#-----------------------------------------------------------------------
# Special internal functions for building tests.
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# function geant4_add_test(<name> COMMAND cmd [arg1... ]
#                          [PRECMD cmd [arg1...]] [POSTCMD cmd [arg1...]]
#                           [OUTPUT outfile] [ERROR errfile]
#                           [WORKING_DIRECTORY directory]
#                           [ENVIRONMENT var1=val1 var2=val2 ...
#                           [DEPENDS test1 ...]
#                           [TIMEOUT seconds]
#                           [DEBUG]
#                           [SOURCE_DIR dir] [BINARY_DIR dir]
#                           [BUILD target] [PROJECT project]
#                           [PASSREGEX exp] [FAILREGEX epx]
#                           [LABELS label1 label2 ...])
#
function(geant4_add_test test)
  if(NOT CMAKE_PROJECT_NAME STREQUAL Geant4)
    message(WARNING "geant4_add_test is only for internal Geant4 usage, test '${test}' is disabled")
    return()
  endif()

  cmake_parse_arguments(ARG
    "DEBUG"
    "TIMEOUT;BUILD;OUTPUT;ERROR;SOURCE_DIR;BINARY_DIR;PROJECT;PASSREGEX;FAILREGEX;WORKING_DIRECTORY"
    "COMMAND;PRECMD;POSTCMD;ENVIRONMENT;DEPENDS;LABELS"
    ${ARGN})

  if(CMAKE_CONFIGURATION_TYPES)
    set(_cfg $<CONFIGURATION>/)
  endif()

  # COMMAND AND BUILD: split test
  # - In this case, we have to create a -build and a -run test with the latter depending on the former
  # NOT COMMAND AND BUILD: pure build
  # COMMAND AND NOT BUILD: pure test
  if(ARG_COMMAND AND ARG_BUILD)
    set(_is_split_test TRUE)
  endif()

  #- Handle COMMAND argument
  list(LENGTH ARG_COMMAND _len)
  if(_len LESS 1)
    if(NOT ARG_BUILD)
      message(FATAL_ERROR "GEANT4_ADD_TEST: command is mandatory (without build)")
    endif()
  else()
    list(GET ARG_COMMAND 0 _prg)
    list(REMOVE_AT ARG_COMMAND 0)
    if(NOT IS_ABSOLUTE ${_prg})
      if(TARGET ${_prg})
        set(_prg "$<TARGET_FILE:${_prg}>")
      else()
        set(_prg ${CMAKE_CURRENT_BINARY_DIR}/${_cfg}${_prg})
      endif()
    elseif(EXISTS ${_prg})
      # Calling a prexisting/system program
    else()
      get_filename_component(_path ${_prg} PATH)
      get_filename_component(_file ${_prg} NAME)
      set(_prg ${_path}/${_cfg}${_file})
    endif()
    set(_cmd ${_prg} ${ARG_COMMAND})
    string(REPLACE ";" "#" _cmd "${_cmd}")
    string(REPLACE "=" "@" _cmd "${_cmd}")
  endif()

  set(_command ${CMAKE_COMMAND} -DTST=${test} -DCMD=${_cmd})

  #- Handle PRE and POST commands
  if(ARG_PRECMD)
    set(_pre ${ARG_PRECMD})
    string(REPLACE ";" "#" _pre "${_pre}")
    set(_command ${_command} -DPRE=${_pre})
  endif()
  if(ARG_POSTCMD)
    set(_post ${ARG_POSTCMD})
    string(REPLACE ";" "#" _post "${_post}")
    set(_command ${_command} -DPOST=${_post})
  endif()

  #- Handle OUTPUT, ERROR, DEBUG arguments
  if(ARG_OUTPUT)
    set(_command ${_command} -DOUT=${ARG_OUTPUT})
  endif()

  if(ARG_ERROR)
    set(_command ${_command} -DERR=${ARG_ERROR})
  endif()

  if(ARG_DEBUG)
    set(_command ${_command} -DDBG=ON)
  endif()

  if(ARG_WORKING_DIRECTORY)
    set(_command ${_command} -DCWD=${ARG_WORKING_DIRECTORY})
  endif()

  if(ARG_TIMEOUT)
    set(_command ${_command} -DTIM=${ARG_TIMEOUT})
  endif()

  #- Handle ENVIRONMENT argument
  if(ARG_ENVIRONMENT)
    string(REPLACE ";" "#" _env "${ARG_ENVIRONMENT}")
    string(REPLACE "=" "@" _env "${_env}")
    set(_command ${_command} -DENV=${_env})
  endif()

  #- Locate the test driver
  # This is one of the few cases where we require ..._SOURCE_DIR to be scoped to Geant4
  # because we may have called geant4_add_test after a project() call, as in the integration tests
  set(_driver ${Geant4_SOURCE_DIR}/cmake/Modules/G4TestDriver.cmake)
  if(NOT EXISTS ${_driver})
    message(FATAL_ERROR "geant4_add_test: G4TestDriver.cmake not found!")
  endif()
  set(_command ${_command} -P ${_driver})

  #- Now we can actually add the test
  if(ARG_BUILD)
    if(NOT ARG_SOURCE_DIR)
      set(ARG_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
    if(NOT ARG_BINARY_DIR)
      set(ARG_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
    endif()
    if(NOT ARG_PROJECT)
       if(NOT PROJECT_NAME STREQUAL "Geant4")
         set(ARG_PROJECT ${PROJECT_NAME})
       else()
         set(ARG_PROJECT ${ARG_BUILD})
       endif()
    endif()

    set(__build_test_name "${test}")
    set(__run_test_name "")

    if(_is_split_test)
      set(__build_test_name "${test}-build")
      set(__run_test_name "${test}")
    endif()

    # Build part of the test
    # Again, note that we scope Geant4_DIR to Geant4_BINARY_DIR so we don't accidentally pickup
    # local or higher scopes
    add_test(NAME ${__build_test_name} COMMAND ${CMAKE_CTEST_COMMAND}
        --build-and-test  ${ARG_SOURCE_DIR} ${ARG_BINARY_DIR}
        --build-generator ${CMAKE_GENERATOR}
        --build-makeprogram ${CMAKE_MAKE_PROGRAM}
        --build-target ${ARG_BUILD}
        --build-project ${ARG_PROJECT}
        --build-config $<CONFIGURATION>
        --build-noclean
        --build-options
          --no-warn-unused-cli
          -DGeant4_DIR=${Geant4_BINARY_DIR}
          -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
          -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
          -DCMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG}
          -DCMAKE_C_FLAGS_MINSIZEREL=${CMAKE_C_FLAGS_MINSIZEREL}
          -DCMAKE_C_FLAGS_RELEASE=${CMAKE_C_FLAGS_RELEASE}
          -DCMAKE_C_FLAGS_RELWITHDEBINFO=${CMAKE_C_FLAGS_RELWITHDEBINFO}
          -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
          -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
          -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
          -DCMAKE_CXX_FLAGS_MINSIZEREL=${CMAKE_CXX_FLAGS_MINSIZEREL}
          -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
          -DCMAKE_CXX_FLAGS_RELWITHDEBINFO=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
          -DCMAKE_DISABLE_FIND_PACKAGE_ROOT=$<BOOL:${CMAKE_DISABLE_FIND_PACKAGE_ROOT}>
    )
    set_property(TEST ${__build_test_name} PROPERTY ENVIRONMENT Geant4_DIR=${Geant4_BINARY_DIR})

    # Build part of the test should have additional regex, and *must* have same labels
    if(ARG_FAILREGEX)
      set_property(TEST ${__build_test_name} PROPERTY FAIL_REGULAR_EXPRESSION "warning:|(${ARG_FAILREGEX})")
    else()
      set_property(TEST ${__build_test_name} PROPERTY FAIL_REGULAR_EXPRESSION "warning:")
    endif()

    if(ARG_LABELS)
      set_property(TEST ${__build_test_name} PROPERTY LABELS ${ARG_LABELS})
    else()
      set_property(TEST ${__build_test_name} PROPERTY LABELS Nightly)
    endif()

    # (Optional) Run part of the test
    if(__run_test_name)
      add_test(NAME ${__run_test_name} COMMAND ${_command})
      set_property(TEST ${__run_test_name} PROPERTY DEPENDS ${__build_test_name})
      if(ARG_FAILREGEX)
        set_property(TEST ${__run_test_name} PROPERTY FAIL_REGULAR_EXPRESSION ${ARG_FAILREGEX})
      endif()

      # If WORKING_DIRECTORY supplied, make sure testdriver is run in
      # that directory, otherwise use specific binary dir
      if(ARG_WORKING_DIRECTORY)
        set_property(TEST ${__run_test_name} PROPERTY WORKING_DIRECTORY "${ARG_WORKING_DIRECTORY}")
      else()
        set_property(TEST ${__run_test_name} PROPERTY WORKING_DIRECTORY "${ARG_BINARY_DIR}")
      endif()
    endif()
  else()
    add_test(NAME ${test} COMMAND ${_command})
    if(ARG_FAILREGEX)
      set_property(TEST ${test} PROPERTY FAIL_REGULAR_EXPRESSION ${ARG_FAILREGEX})
    endif()

    # If WORKING_DIRECTORY supplied, make sure testdriver is run in
    # that directory.
    if(ARG_WORKING_DIRECTORY)
      set_property(TEST ${test} PROPERTY WORKING_DIRECTORY "${ARG_WORKING_DIRECTORY}")
    endif()
  endif()

  #- Handle TIMOUT and DEPENDS arguments
  if(ARG_TIMEOUT)
    set_property(TEST ${test} PROPERTY TIMEOUT ${ARG_TIMEOUT})
  endif()

  if(ARG_DEPENDS)
    # MUST be append because a split test may set this earlier
    set_property(TEST ${test} APPEND PROPERTY DEPENDS ${ARG_DEPENDS})
  endif()

  if(ARG_PASSREGEX)
    set_property(TEST ${test} PROPERTY PASS_REGULAR_EXPRESSION ${ARG_PASSREGEX})
  endif()

  if(ARG_LABELS)
    set_property(TEST ${test} PROPERTY LABELS ${ARG_LABELS})
  else()
    set_property(TEST ${test} PROPERTY LABELS Nightly)
  endif()
endfunction()

