#-----------------------------------------------------------------------
# Geant4 test driver
#   Script arguments: 
#     CMD command to be executed for the test
#     PRE command to be executed before the test command
#     POST command to be executed after the test command
#     OUT file to collect stdout
#     ERR file to collect stderr
#     ENV evironment VAR1=Value1;VAR2=Value2
#     CWD current working directory
#     TST test name (used to name output/error files)
#     TIM timeout 
#     DBG debug flag

if(DBG)
  message(STATUS "ENV=${ENV}")
endif()

#-----------------------------------------------------------------------
# Message arguments
#
if(CMD)
  string(REPLACE "#" ";" _cmd "${CMD}")
  string(REPLACE "@" "=" _cmd "${_cmd}")
  if(DBG)
    set(_cmd gdb --args ${_cmd})
  endif()
endif()

if(PRE)
  string(REPLACE "#" ";" _pre ${PRE})
endif()

if(POST)
  string(REPLACE "#" ";" _post ${POST})
endif()

if(OUT)
  set(_out OUTPUT_FILE ${OUT})
else()
  if(NOT DBG)
    set(_out OUTPUT_VARIABLE _outvar)
  endif()
endif()

if(ERR)
  set(_err ERROR_FILE ${ERR})
else()
  if(NOT DBG)
    set(_err ERROR_VARIABLE _errvar)
  endif()
endif()

if(TIM)
  math(EXPR _timeout "${TIM} - 120")
else()
  if(NOT $ENV{CTEST_TIMEOUT} STREQUAL "")
     math(EXPR _timeout "$ENV{CTEST_TIMEOUT} - 120")
  else()	  
     math(EXPR _timeout "1500 - 120")
  endif()	  
endif()

if(CWD)
  set(_cwd WORKING_DIRECTORY ${CWD})
endif()

#-----------------------------------------------------------------------
# Environment settings
#
if(ENV)
  string(REPLACE "@" "=" _env ${ENV})
  string(REPLACE "#" ";" _env ${_env})
  foreach(pair ${_env})
    string(REPLACE "=" ";" pair ${pair})
    list(GET pair 0 var)
    list(GET pair 1 val)
    set(ENV{${var}} ${val})
    if(DBG)
      message(STATUS "testdriver[ENV]:${var}==>${val}")
    endif()
  endforeach()
endif()

#-----------------------------------------------------------------------
# Execute pre command
#
if(PRE)
  execute_process(COMMAND ${_pre} ${_cwd} RESULT_VARIABLE _rc)
  if(_rc)
    message(FATAL_ERROR "pre-command error code : ${_rc}")
  endif()
endif()

#-----------------------------------------------------------------------
# Execute test
#
if(CMD)
  execute_process(COMMAND ${_cmd} ${_out} ${_err} ${_cwd} TIMEOUT ${_timeout} RESULT_VARIABLE _rc)
  message("G4Test rc: ${_rc}")
  if(_errvar)
    message("G4Test stderr:\n ${_errvar}")
    if(TST)
      file(WRITE ${TST}.err ${_errvar})
    endif()
  endif()
  if(_outvar)
    message("G4Test stdout:\n ${_outvar}")
    if(TST)
      file(WRITE ${TST}.out ${_outvar})
    endif()
  endif()

  # - FATAL_ERROR if the test returned an error code or wrote anything to
  # the stderr
  if(_errvar OR _rc)
    message(FATAL_ERROR "Test failed!!!")
  endif()
endif()

#-----------------------------------------------------------------------
# Execute post test command
#
if(POST)
  execute_process(COMMAND ${_post} ${_cwd} RESULT_VARIABLE _rc)
  if(_rc)
    message(FATAL_ERROR "post-command error code : ${_rc}")
  endif()
endif()

