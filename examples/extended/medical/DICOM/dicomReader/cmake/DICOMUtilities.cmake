# Utility macros for DICOM example
# called from the DICOM example CMakeLists.txt
#

# - Include guard
if(__dicomutilities_isloaded)
  return()
endif()
set(__dicomutilities_isloaded YES)

# - for cmake_parse_arguments
include(CMakeParseArguments)

macro(DICOM_BUILD_LIBRARY)

    set(_options    )                       # options
    set(_onevalue   BUILD_SHARED            # single-value
                    BUILD_STATIC
                    OUTPUT_NAME
                    TARGET_NAME)
    set(_multival   SOURCES                 # multi-value
                    LINK_LIBRARIES
                    COMPILE_DEFINITIONS)

    cmake_parse_arguments(
        LIBRARY "${_options}" "${_onevalue}" "${_multival}" ${ARGN})

    # static library
    if(LIBRARY_BUILD_STATIC)
        # add static library
        add_library(${LIBRARY_TARGET_NAME}-static STATIC ${LIBRARY_SOURCES})
        # link
        target_link_libraries(${LIBRARY_TARGET_NAME}-static ${LIBRARY_LINK_LIBRARIES})
        # properties
        set_target_properties(${LIBRARY_TARGET_NAME}-static PROPERTIES
            OUTPUT_NAME ${LIBRARY_OUTPUT_NAME}
            POSITION_INDEPENDENT_CODE OFF)
        # append this to list of libraries to install
        list(APPEND ${PROJECT_NAME}_INSTALL_LIBRARIES ${LIBRARY_TARGET_NAME}-static)
    endif(LIBRARY_BUILD_STATIC)

    # shared library
    if(LIBRARY_BUILD_SHARED)
        ####################################################
        #   This handles Windows issues
        set(_TYPE SHARED )
        set(_EXTRA_PROPS )
        # CMake >= 3.4 has WINDOWS_EXPORT_ALL_SYMBOLS
        if(WIN32 AND (CMAKE_VERSION VERSION_LESS 3.4))
            # don't bother with dll
            set(_TYPE STATIC)
        elseif(WIN32)
            set(_EXTRA_PROPS WINDOWS_EXPORT_ALL_SYMBOLS ON)
        endif()
        ####################################################
        # add shared (or static if Windows and CMake < 3.4)
        add_library(${LIBRARY_TARGET_NAME} ${_TYPE} ${LIBRARY_SOURCES})
        # link
        target_link_libraries(${LIBRARY_TARGET_NAME} ${LIBRARY_LINK_LIBRARIES})
        # properties
        set_target_properties(${LIBRARY_TARGET_NAME} PROPERTIES
            OUTPUT_NAME ${LIBRARY_OUTPUT_NAME}
            POSITION_INDEPENDENT_CODE ON
            ${_EXTRA_PROPS})
        # append this to list of libraries to install
        list(APPEND ${PROJECT_NAME}_INSTALL_LIBRARIES ${LIBRARY_TARGET_NAME})
        # cleanup
        unset(_TYPE)
        unset(_EXTRA_PROPS)
    endif(LIBRARY_BUILD_SHARED)

    add_library(${PROJECT_NAME}::library ALIAS ${LIBRARY_TARGET_NAME}${_geant4_lib_use_suffix})

endmacro(DICOM_BUILD_LIBRARY)
