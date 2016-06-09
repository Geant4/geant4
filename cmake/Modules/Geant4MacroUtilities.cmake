# Geant4MacroUtilities - useful macros and functions for generic tasks
#
# CMake Extensions
# ----------------
# macro set_ifnot(<var> <value>)
#       If variable var is not set, set its value to that provided
#
# function enum_option(<option>
#                      VALUES <value1> ... <valueN>
#                      TYPE   <valuetype>
#                      DOC    <docstring>
#                      [DEFAULT <elem>]
#                      [CASE_INSENSITIVE]
#          Declare a cache variable <option> that can only take values
#          listed in VALUES. TYPE may be FILEPATH, PATH or STRING.
#          <docstring> should describe that option, and will appear in
#          the interactive CMake interfaces. If DEFAULT is provided,
#          <elem> will be taken as the zero-indexed element in VALUES
#          to which the value of <option> should default to if not
#          provided. Otherwise, the default is taken as the first
#          entry in VALUES. If CASE_INSENSITIVE is present, then
#          checks of the value of <option> against the allowed values
#          will ignore the case when performing string comparison.
#
#
# General Geant4
# --------------
# macro geant4_add_feature(<NAME> <DOCSTRING>)
#       Add a Geant4 feature, whose activation is specified by the
#       existence of the variable <NAME>, to the list of enabled/disabled
#       features, plus a docstring describing the feature
#
# macro geant4_print_enabled_features()
#       Print enabled Geant4 features plus their docstrings.
#
# Datasets
# --------
# TODO: Move to dedicated datasets module
# function geant4_latest_version(<dir> <name> <output variable>)
#          Locate latest version of dataset <name> in <dir>, setting value
#          of output variable to the full path to the dataset
#
# Testing
# -------
# TODO: Move to dedicated tests module
# function geant4_add_unit_tests(test1 test2 ... [dir1 ...]
#                                INCLUDE_DIRS dir1 dir2 ...
#                                LIBRARIES library1 library2 ...)
#

include(Geant4MacroDefineModule)
include(Geant4MacroLibraryTargets)

#-----------------------------------------------------------------------
# CMAKE EXTENSIONS
#-----------------------------------------------------------------------
# macro set_ifnot(<var> <value>)
#       If variable var is not set, set its value to that provided
#
macro(set_ifnot _var _value)
  if(NOT ${_var})
    set(${_var} ${_value})
  endif()
endmacro()

#-----------------------------------------------------------------------
# function enum_option(<option>
#                      VALUES <value1> ... <valueN>
#                      TYPE   <valuetype>
#                      DOC    <docstring>
#                      [DEFAULT <elem>]
#                      [CASE_INSENSITIVE]
#          Declare a cache variable <option> that can only take values
#          listed in VALUES. TYPE may be FILEPATH, PATH or STRING.
#          <docstring> should describe that option, and will appear in
#          the interactive CMake interfaces. If DEFAULT is provided,
#          <elem> will be taken as the zero-indexed element in VALUES
#          to which the value of <option> should default to if not
#          provided. Otherwise, the default is taken as the first
#          entry in VALUES. If CASE_INSENSITIVE is present, then
#          checks of the value of <option> against the allowed values
#          will ignore the case when performing string comparison.
#
function(enum_option _var)
  set(options CASE_INSENSITIVE)
  set(oneValueArgs DOC TYPE DEFAULT)
  set(multiValueArgs VALUES)
  cmake_parse_arguments(_ENUMOP "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # - Validation as needed arguments 
  if(NOT _ENUMOP_VALUES)
    message(FATAL_ERROR "enum_option must be called with non-empty VALUES\n(Called for enum_option '${_var}')")
  endif()

  # - Set argument defaults as needed
  if(_ENUMOP_CASE_INSENSITIVE)
    set(_ci_values )
    foreach(_elem ${_ENUMOP_VALUES})
      string(TOLOWER "${_elem}" _ci_elem)
      list(APPEND _ci_values "${_ci_elem}")
    endforeach()
    set(_ENUMOP_VALUES ${_ci_values})
  endif()

  set_ifnot(_ENUMOP_TYPE STRING)
  set_ifnot(_ENUMOP_DEFAULT 0)
  list(GET _ENUMOP_VALUES ${_ENUMOP_DEFAULT} _default)

  if(NOT DEFINED ${_var})
    set(${_var} ${_default} CACHE ${_ENUMOP_TYPE} "${_ENUMOP_DOC} (${_ENUMOP_VALUES})")
  else()
    set(_var_tmp ${${_var}})
    if(_ENUMOP_CASE_INSENSITIVE)
      string(TOLOWER ${_var_tmp} _var_tmp)
    endif() 

    list(FIND _ENUMOP_VALUES ${_var_tmp} _elem)
    if(_elem LESS 0)
      message(FATAL_ERROR "Value '${${_var}}' for variable ${_var} is not allowed\nIt must be selected from the set: ${_ENUMOP_VALUES} (DEFAULT: ${_default})\n")
    else()
      # - convert to lowercase
      if(_ENUMOP_CASE_INSENSITIVE)
        set(${_var} ${_var_tmp} CACHE ${_ENUMOP_TYPE} "${_ENUMOP_DOC} (${_ENUMOP_VALUES})" FORCE)
      endif()
    endif()
  endif()
endfunction()

#-----------------------------------------------------------------------
# GENERAL GEANT4
#-----------------------------------------------------------------------
# macro geant4_add_feature(<NAME> <DOCSTRING>)
#       Add a Geant4 feature, whose activation is specified by the
#       existence of the variable <NAME>, to the list of enabled/disabled
#       features, plus a docstring describing the feature
#
macro(GEANT4_ADD_FEATURE _var _description)
  if(${_var})
    set_property(GLOBAL APPEND PROPERTY GEANT4_ENABLED_FEATURES ${_var})
  else(${_var})
    set_property(GLOBAL APPEND PROPERTY GEANT4_DISABLED_FEATURES ${_var})
  endif(${_var})

  set_property(GLOBAL PROPERTY ${_var}_DESCRIPTION "${_description}")
endmacro()

#-----------------------------------------------------------------------
# macro geant4_print_enabled_features()
#       Print enabled Geant4 features plus their docstrings.
#
macro(GEANT4_PRINT_ENABLED_FEATURES)
  set(_currentFeatureText "The following Geant4 features are enabled:")
  get_property(_enabledFeatures GLOBAL PROPERTY GEANT4_ENABLED_FEATURES)

  foreach(_feature ${_enabledFeatures})
    set(_currentFeatureText "${_currentFeatureText}\n${_feature}")

    get_property(_desc GLOBAL PROPERTY ${_feature}_DESCRIPTION)

    if(_desc)
      set(_currentFeatureText "${_currentFeatureText}: ${_desc}")
      set(_desc NOTFOUND)
    endif(_desc)
  endforeach(_feature)

  message(STATUS "${_currentFeatureText}\n")
endmacro()

#-----------------------------------------------------------------------
# GEANT4 DATASETS
#-----------------------------------------------------------------------
# function geant4_latest_version(<dir> <name> <output variable>)
#          Locate latest version of dataset <name> in <dir>, setting value
#          of output variable to the full path to the dataset
#
function(GEANT4_LATEST_VERSION dir name var)
  file(GLOB files RELATIVE ${dir} ${dir}/${name}*)
  set(newer)
  foreach(file ${files})
    string(REPLACE ${name} "" version ${file})
    if("${version}" VERSION_GREATER "${newer}")
      set(newer ${version})
    endif()
  endforeach()
  set(${var} ${dir}/${name}${newer} PARENT_SCOPE)
endfunction()


#-----------------------------------------------------------------------
# GEANT4 TESTING
#-----------------------------------------------------------------------
# function geant4_add_unit_tests(test1 test2 ... [dir1 ...]
#                                INCLUDE_DIRS dir1 dir2 ...
#                                LIBRARIES library1 library2 ...)
#
function(GEANT4_ADD_UNIT_TESTS)
  CMAKE_PARSE_ARGUMENTS(ARG "" "" "INCLUDE_DIRS;LIBRARIES" ${ARGN})

  foreach(incdir ${ARG_INCLUDE_DIRS})
    if(IS_ABSOLUTE ${incdir})
      include_directories(${incdir})
    else()
      include_directories(${CMAKE_SOURCE_DIR}/source/${incdir})
    endif()
  endforeach()
  
  if(ARG_UNPARSED_ARGUMENTS)
    set(tnames ${ARG_UNPARSED_ARGUMENTS})
  else()
    set(tnames test*.cc)
  endif()
  
  set(alltests)
  foreach(tname ${tnames})
    if(tname STREQUAL ".")
      set(tests ".")
    else()
      file(GLOB tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${tname})
    endif()
    set(alltests ${alltests} ${tests})
  endforeach()
  
  if(NOT TARGET tests)
    add_custom_target(tests)
  endif()
  
  foreach(test ${alltests})
    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${test})
      file(GLOB sources ${test}/src/*.cc)
      include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${test}/include)
      file(GLOB test ${test}/*.cc)
    else()
      set(sources)
    endif()
    get_filename_component(name ${test} NAME_WE)
    add_executable(${name} EXCLUDE_FROM_ALL ${test} ${sources})
    target_link_libraries(${name} ${ARG_LIBRARIES})
    set_target_properties(${name} PROPERTIES OUTPUT_NAME ${name})
    add_dependencies(tests ${name})
    add_test(NAME ${name} COMMAND ${name})
    set_property(TEST ${name} PROPERTY LABELS UnitTests)
    set_property(TEST ${name} PROPERTY TIMEOUT 60)
  endforeach()
endfunction()
