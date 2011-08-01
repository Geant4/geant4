# - Use file for Geant4
# This file should be included after a find_package call has successfully
# located Geant4.
#
#  Geant4_USE_FILE : Point to the location of the use file for the found
#                    Geant4 installation.
#
# Inclusion of this file results in the addition of the Geant4 compile
# definitions and include directories to those of the directory in which
# this files is included.
#
# The recommended Geant4 compiler flags are also added to CMAKE_CXX_FLAGS
# but duplicated flags are NOT removed.
# Advanced users requiring special sets of flags, or the removal of duplicate
# flags should therefore *not* use this file, preferring the direct use of the 
# Geant4_XXXX variables set by the Geant4Config file.
#
# The use file also defines a simple macro to help with collating sources
# for simple user applications.
#
#  macro GEANT4_COLLATE_APPLICATION_SOURCES(sources)
#        Create the list of sources in a standard Geant4 application example.
#        If the project is organized as
#        +- project
#           +- CMakeLists.txt
#           +- project.cc
#           +- include/
#           |  +- project_impl_a.hh
#           |  +- ...
#           +- src/
#              +- project_impl_a.cc
#              +- ...
#        Then when called in the CMakeLists.txt under 'project', it will
#        collate all the .hh and .cc files into the sources variable, and
#        add the include directory to those to be searched for.
#        The CMakeLists.txt for project would then read:
#       
#        cmake_minimum_required(VERSION 2.6.4)
#        project(myproject)
#
#        find_package(Geant4 REQUIRED)
#        include(${Geant4_USE_FILE})
#
#        GEANT4_COLLATE_APPLICATION_SOURCES(my_project_sources)
#
#        add_executable(myproject ${my_project_sources})
#        target_link_libraries(myproject ${Geant4_LIBRARIES})
#
#
#        This gives you maximum flexibility in setting up your project, as
#        you can either use this macro for simplicity and to easily convert
#        an existing application to CMake, or you can set the sources manually
#        with your own arrangement of headers and source files.
#

#----------------------------------------------------------------------------
# macro GEANT4_COLLATE_SOURCES(source_dest_var)
#
macro(GEANT4_COLLATE_APPLICATION_SOURCES source_dest_var)
    file(GLOB_RECURSE 
        ${source_dest_var} 
        ${CMAKE_CURRENT_SOURCE_DIR}/*.hh 
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cc
    )
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
endmacro()


#----------------------------------------------------------------------------
# We need to set the compile definitions and include directories
#
add_definitions(${Geant4_DEFINITIONS})
include_directories(${Geant4_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Because Geant4 is sensitive to the compiler flags, lets set the base set
# here. This reproduces as far as possible the behaviour of the original
# makefile system.
#
set(CMAKE_CXX_FLAGS                "${Geant4_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "${Geant4_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${Geant4_CXX_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE        "${Geant4_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${Geant4_CXX_FLAGS_RELWITHDEBINFO}")



#----------------------------------------------------------------------------
# Special internal functions for building tests.
#
include(CMakeMacroParseArguments)

#----------------------------------------------------------------------------
# function GEANT4_LINK_LIBRARY( <name> source1 source2 ...
#                               [TYPE STATIC|SHARED] 
#                               LIBRARIES library1 library2 ... )
#
function(GEANT4_LINK_LIBRARY library)
  CMAKE_PARSE_ARGUMENTS(ARG "TYPE;LIBRARIES" "" ${ARGN})
  set(sources)

  # - Fill sources
  foreach(fp ${ARG_UNPARSED_ARGUMENTS})  
    if(IS_ABSOLUTE ${fp}) 
      file(GLOB files ${fp})     
    else()
      file(GLOB files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${fp})
    endif()
    if(files) 
      set(sources ${sources} ${files})
    else()
      set(sources ${sources} ${fp})
    endif()
  endforeach()

  # - Shared library unless specified
  if(NOT ARG_TYPE)
    set(ARG_TYPE SHARED)
  endif()

  # - Make sure we can access our own headers
  include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include)

  # - Deal with Win32 DLLs that don't export via declspec
  if(WIN32 AND ARG_TYPE STREQUAL SHARED)
    # - Dummy archive library
    add_library( ${library}-arc STATIC EXCLUDE_FROM_ALL ${sources})

    # - Use genwindef to create .def file listing symbols
    add_custom_command(
      OUTPUT ${library}.def
      COMMAND ${genwindef_cmd} -o ${library}.def -l ${library} ${LIBRARY_OUTPUT_PATH}/${CMAKE_CFG_INTDIR}/${library}-arc.lib
      DEPENDS ${library}-arc genwindef)

    #- Dummy cpp file needed to satisfy Visual Studio.
    file( WRITE ${CMAKE_CURRENT_BINARY_DIR}/${library}.cpp "// empty file\n" )
    add_library( ${library} SHARED ${library}.cpp ${library}.def)
    target_link_libraries(${library} ${library}-arc ${ARG_LIBRARIES})
    set_target_properties(${library} PROPERTIES LINK_INTERFACE_LIBRARIES ${ARG_LIBRARIES} ${Geant4_LIBRARIES})
  else()
    add_library( ${library} ${ARG_TYPE} ${sources})
    target_link_libraries(${library} ${ARG_LIBRARIES} ${Geant4_LIBRARIES})
  endif()

endfunction()


#----------------------------------------------------------------------------
# function GEANT4_EXECUTABLE( <name> source1 source2 ...
#                             LIBRARIES library1 library2 ... )
#
function(GEANT4_EXECUTABLE executable)
  CMAKE_PARSE_ARGUMENTS(ARG "" "" "LIBRARIES" ${ARGN})
  set(sources)
  foreach( fp ${ARG_UNPARSED_ARGUMENTS})
    file(GLOB files ${fp})
    if(files)
      set( sources ${sources} ${files})
    else()
      set( sources ${sources} ${fp})
    endif()
  endforeach()
  include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include ${GEANT4_INCLUDE_DIR})
  add_executable(${executable} EXCLUDE_FROM_ALL ${sources})
  target_link_libraries(${executable} ${ARG_LIBRARIES} ${Geant4_LIBRARIES} )
  set_target_properties(${executable} PROPERTIES OUTPUT_NAME ${executable})

endfunction()


#----------------------------------------------------------------------------
# function GEANT4_ADD_TEST( <name> COMMAND cmd [arg1... ] 
#                           [OUTPUT outfile] [ERROR errfile]
#                           [ENVIRONMENT var1=val1 var2=val2 ...
#                           [TIMEOUT seconds] 
#                           [DEBUG]
#                           [BUILD target] )
#
function(GEANT4_ADD_TEST test)
  if(NOT CMAKE_PROJECT_NAME STREQUAL Geant4)
    message(WARNING "function GEANT4_ADD_TEST is only for internal Geant4 usage")
    return()
  endif()
  CMAKE_PARSE_ARGUMENTS(ARG "DEBUG" "TIMEOUT;BUILD;OUTPUT;ERROR" "COMMAND;ENVIRONMENT;DEPENDS" ${ARGN})

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
      set(_prg ${CMAKE_CURRENT_BINARY_DIR}/${_prg})
    endif()
    set(_cmd ${_prg} ${ARG_COMMAND})
    string(REPLACE ";" ":" _cmd "${_cmd}")
  endif()

  set(_command ${CMAKE_COMMAND} -DCMD=${_cmd})

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

  #- Handle ENVIRONMENT argument
  if(ARG_ENVIRONMENT)
    string(REPLACE ";" ":" _env "${ARG_ENVIRONMENT}")
    string(REPLACE "=" "@" _env "${_env}")
    set(_command ${_command} -DENV=${_env})
  endif()

  #- Locate the test driver
  set(_driver ${CMAKE_SOURCE_DIR}/cmake/Modules/Geant4TestDriver.cmake)
  if(NOT EXISTS ${_driver})
    message(FATAL_ERROR "GEANT4_ADD_TEST: Geant4TestDriver.cmake not found!")
  endif()
  set(_command ${_command} -P ${_driver})


  #- Now we can actually add the test
  if(ARG_BUILD)
    add_test(NAME ${test} COMMAND ${CMAKE_CTEST_COMMAND}
      --build-and-test  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${ARG_BUILD}
      --build-noclean
      --test-command ${_command} )
    set_property(TEST ${test} PROPERTY ENVIRONMENT Geant4_DIR=${CMAKE_BINARY_DIR})
  else()
    add_test(NAME ${test} COMMAND ${_command})
  endif()

  #- Handle TIMOUT and DEPENDS arguments
  if(ARG_TIMEOUT)
    set_property(TEST ${test} PROPERTY TIMEOUT ${ARG_TIMEOUT})
  endif()

  if(ARG_DEPENDS)
    set_property(TEST ${test} PROPERTY DEPENDS ${ARG_DEPENDS})
  endif()

endfunction()

