# PTLCMakeUtilities - useful macros and functions for generic tasks
#
# CMake Extensions
# ----------------

# * Include guard
include_guard(DIRECTORY)
include(CMakeParseArguments)

# -------------------------------------------------------------------------------------- #
# function ptl_message_on_change(<VAR> <MESSAGE>) Watch <VAR> and print <MESSAGE> and
# <VAR>'s value when the value changes
#
function(ptl_message_on_change _var _prefix)
    set(_ptl_watch_var_name PTL_WATCH_VALUE_${_var})
    if(DEFINED ${_ptl_watch_var_name})
        if("${${_var}}" STREQUAL "${${_ptl_watch_var_name}}")
            return()
        endif()
    endif()

    set(${_ptl_watch_var_name}
        "${${_var}}"
        CACHE INTERNAL "Last value of ${_var}" FORCE)
    message(STATUS "${_prefix} - ${${_var}}")
endfunction()

# -------------------------------------------------------------------------------------- #
# function ptl_add_feature(<NAME> <DOCSTRING>) Store activation of a feature and a
# docstring describing it. <NAME> should be a variable whose value indicates activation of
# the feature
#
function(ptl_add_feature _var _description)
    if(${_var})
        set_property(GLOBAL APPEND PROPERTY PTL_ENABLED_FEATURES ${_var})
    endif()

    set_property(GLOBAL PROPERTY PTL_FEATURE_${_var}_DESCRIPTION "${_description}")
endfunction()

# -----------------------------------------------------------------------
# function ptl_print_features() Print enabled features plus their docstrings.
#
function(ptl_print_features)
    set(_currentFeatureText "The following PTL features are enabled:")
    get_property(_enabledFeatures GLOBAL PROPERTY PTL_ENABLED_FEATURES)

    foreach(_feature ${_enabledFeatures})
        set(_currentFeatureText "${_currentFeatureText}\n   - ${_feature}")

        get_property(_desc GLOBAL PROPERTY PTL_FEATURE_${_feature}_DESCRIPTION)

        if(_desc)
            set(_currentFeatureText "${_currentFeatureText}: ${_desc}")
            set(_desc NOTFOUND)
        endif()
    endforeach()

    message(STATUS "${_currentFeatureText}\n")
endfunction()

# -------------------------------------------------------------------------------------- #
# function ptl_add_option(<OPTION_NAME> <DOCSTRING> <DEFAULT_SETTING>) Add an option for
# master project only, always reporting its value
#
function(ptl_add_option _NAME _MESSAGE _DEFAULT)
    if(PTL_MASTER_PROJECT)
        option(${_NAME} "${_MESSAGE}" ${_DEFAULT})
        ptl_add_feature(${_NAME} "${_MESSAGE}")
    endif()
    # Always on configure time reporting
    ptl_message_on_change(${_NAME} "Building PTL with option ${_NAME}")
endfunction()

# ----------------------------------------------------------------------------
# function ptl_build_library Encapsulates common build/install for main PTL shared/static
# libs
#
function(ptl_build_library)
    cmake_parse_arguments(LIB "" "TYPE;TARGET_NAME;OUTPUT_NAME" "SOURCES;EXTRA_ARGS"
                          ${ARGN})

    if(NOT LIB_OUTPUT_NAME)
        set(LIB_OUTPUT_NAME ${LIB_TARGET_NAME})
    endif()

    add_library(${LIB_TARGET_NAME} ${LIB_TYPE} ${PTL_EXCLUDE_FROM_ALL})

    add_library(${PROJECT_NAME}::${LIB_TARGET_NAME} ALIAS ${LIB_TARGET_NAME})

    target_sources(${LIB_TARGET_NAME} PRIVATE ${LIB_SOURCES})

    target_compile_definitions(${LIB_TARGET_NAME} PRIVATE $<$<CONFIG:Debug>:DEBUG>)

    target_compile_features(${LIB_TARGET_NAME} PUBLIC cxx_std_${CMAKE_CXX_STANDARD})

    # Subproject overrides only if not a master project
    if(NOT PTL_MASTER_PROJECT)
        target_compile_options(
            ${LIB_TARGET_NAME}
            PRIVATE $<$<COMPILE_LANGUAGE:C>:${${PROJECT_NAME}_C_FLAGS}>
                    $<$<COMPILE_LANGUAGE:CXX>:${${PROJECT_NAME}_CXX_FLAGS}>)
    endif()

    set_target_properties(
        ${LIB_TARGET_NAME}
        PROPERTIES OUTPUT_NAME ${LIB_OUTPUT_NAME}
                   VERSION ${${PROJECT_NAME}_VERSION}
                   SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR}
                   WINDOWS_EXPORT_ALL_SYMBOLS ON)

    # Install the targets and export libraries
    if(NOT "${LIB_TYPE}" STREQUAL "OBJECT")
        install(
            TARGETS ${LIB_TARGET_NAME}
            EXPORT ${PROJECT_NAME}Targets
            COMPONENT Development
            ARCHIVE DESTINATION ${PTL_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${PTL_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${PTL_INSTALL_BINDIR} OPTIONAL)
    endif()
endfunction()
