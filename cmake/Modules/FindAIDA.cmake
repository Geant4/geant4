# - Finds AIDA instalation
# This module sets up AIDA information 
# It defines:
# AIDA_FOUND           If the AIDA is found
# AIDA_INCLUDE_DIR     PATH to the include directory
# AIDA_LIBRARIES       AIDA installation libraries
# AIDA_IMPLEMENTATION  The implementation of AIDA found

find_program(AIDA_CONFIG_EXECUTABLE aida-config)

if(NOT AIDA_CONFIG_EXECUTABLE)
  set(AIDA_FOUND FALSE)
else()    
  set(AIDA_FOUND TRUE)

  execute_process(
    COMMAND ${AIDA_CONFIG_EXECUTABLE} --include
    OUTPUT_VARIABLE AIDA_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
    #remove -I prefix 
    string(REGEX REPLACE "-I" "" AIDA_INCLUDE_DIR "${AIDA_INCLUDE_DIR}")

  execute_process(
    COMMAND ${AIDA_CONFIG_EXECUTABLE} --libs
    OUTPUT_VARIABLE AIDA_LIBRARIES
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(
    COMMAND ${AIDA_CONFIG_EXECUTABLE} --version 
    OUTPUT_VARIABLE AIDA_VERSION
    ERROR_VARIABLE AIDA_VERSION_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(
    COMMAND ${AIDA_CONFIG_EXECUTABLE} --implementation
    OUTPUT_VARIABLE AIDA_IMPLEMENTATION
    ERROR_VARIABLE AIDA_IMPLEMENTATION_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Make variables changeble to the advanced user
  mark_as_advanced(AIDA_CONFIG_EXECUTABLE)

  if(NOT AIDA_FIND_QUIETLY)
    set (OUTPUT "Found AIDA")
    if (AIDA_VERSION)
      set (OUTPUT "${OUTPUT} ${AIDA_VERSION}")
    endif()  
    if (AIDA_IMPLEMENTATION)
      set (OUTPUT "${OUTPUT} implemented by ${AIDA_IMPLEMENTATION}")
    endif()  
    message(STATUS "${OUTPUT}")
  endif()
  
endif()
  
