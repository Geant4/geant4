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

    set(PREFER_SHARED ON)
    if(NOT "${_geant4_lib_use_suffix}" STREQUAL "" AND Geant4_static_FOUND)
        set(PREFER_SHARED OFF)
    endif(NOT "${_geant4_lib_use_suffix}" STREQUAL "" AND Geant4_static_FOUND)

    # static library
    if(LIBRARY_BUILD_STATIC)
        add_library(${LIBRARY_TARGET_NAME}-static STATIC ${LIBRARY_SOURCES})
        target_link_libraries(${LIBRARY_TARGET_NAME}-static ${LIBRARY_LINK_LIBRARIES})
        set_target_properties(${LIBRARY_TARGET_NAME}-static PROPERTIES
            OUTPUT_NAME                 ${LIBRARY_OUTPUT_NAME}
            POSITION_INDEPENDENT_CODE   OFF)
        list(APPEND ${PROJECT_NAME}_INSTALL_LIBRARIES
            ${LIBRARY_TARGET_NAME}-static)
        set(${LIBRARY_TARGET_NAME}-target ${LIBRARY_TARGET_NAME}-static)
    endif(LIBRARY_BUILD_STATIC)

    # shared library
    # do this after static so -target, suffix are preferred when shared built
    # except is PREFER_SHARED is OFF
    if(LIBRARY_BUILD_SHARED)
        add_library(${LIBRARY_TARGET_NAME} SHARED ${LIBRARY_SOURCES})
        target_link_libraries(${LIBRARY_TARGET_NAME} ${LIBRARY_LINK_LIBRARIES})
        set_target_properties(${LIBRARY_TARGET_NAME} PROPERTIES
            OUTPUT_NAME                 ${LIBRARY_OUTPUT_NAME}
            POSITION_INDEPENDENT_CODE   ON)
        list(APPEND ${PROJECT_NAME}_INSTALL_LIBRARIES
            ${LIBRARY_TARGET_NAME})
        if(PREFER_SHARED)
        endif(PREFER_SHARED)
    endif(LIBRARY_BUILD_SHARED)

    set(${LIBRARY_TARGET_NAME}-target ${LIBRARY_TARGET_NAME}${_geant4_lib_use_suffix})
    unset(PREFER_SHARED)

endmacro(DICOM_BUILD_LIBRARY)
