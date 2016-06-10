# - Find ROOT library
# This module sets up ROOT information
# It defines:
# ROOT_FOUND          If the ROOT is found
# ROOT_INCLUDE_DIR    PATH to the include directory
# ROOT_LIBRARIES      Most common libraries
# ROOT_LIBRARY_DIR    PATH to the library directory

find_program(ROOT_CONFIG_EXECUTABLE root-config
             PATHS $ENV{ROOTSYS}/bin)

if(NOT ROOT_CONFIG_EXECUTABLE)
  set(ROOT_FOUND FALSE)
  message(STATUS "NOT Found ROOT.")

else()
  set(ROOT_FOUND TRUE)

  execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --prefix
                  OUTPUT_VARIABLE ROOTSYS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --version
                  OUTPUT_VARIABLE ROOT_VERSION
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --incdir
                  OUTPUT_VARIABLE ROOT_INCLUDE_DIR
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --libs
                  OUTPUT_VARIABLE ROOT_LIBRARIES
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  set(ROOT_LIBRARY_DIR ${ROOTSYS}/lib)

  message(STATUS "Found ROOT ${ROOT_VERSION} in ${ROOTSYS}")

endif()
