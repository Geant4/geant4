# Geant4MacroUtilities - this module defines helper macros and functions
#
# GEANT4_ADD_FEATURE(NAME DESCRIPTION)
#   Use this macro to add a Geant4 specific feature NAME, assumed to be a 
#   boolean, to the list of enabled/disabled features, together with a short 
#   DESCRIPTION.
#
# GEANT4_PRINT_ENABLED_FEATURES()
#   Prints list of enabled Geant4 features and their description only. Just a
#   simplified version of that in FeatureSummary.



macro(GEANT4_ADD_FEATURE _var _description)
    if(${_var})
        set_property(GLOBAL APPEND PROPERTY GEANT4_ENABLED_FEATURES ${_var})
    else(${_var})
        set_property(GLOBAL APPEND PROPERTY GEANT4_DISABLED_FEATURES ${_var})
    endif(${_var})

    set_property(GLOBAL PROPERTY ${_var}_DESCRIPTION "${_description}")
endmacro(GEANT4_ADD_FEATURE)


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
endmacro(GEANT4_PRINT_ENABLED_FEATURES)


#----------------------------------------------------------------------------
#---Helper function to locate the most recent version of a dataset-----------
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


#----------------------------------------------------------------------------
# function GEANT4_ADD_UNIT_TESTS(test1 test2 ... [dir1 ...]
#                                INCLUDE_DIRS dir1 dir2 ...
#                                LIBRARIES library1 library2 ... )
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
