################################################################################
#
#        Creates a 'format' target that runs clang-format
#
################################################################################

find_program(CLANG_FORMATTER
    NAMES
        clang-format-6
        clang-format-6.0
        clang-format-mp-6.0 # macports
        clang-format)

if(CLANG_FORMATTER)
    set(_Source_DIR     ${PROJECT_SOURCE_DIR}/src)
    set(_Header_DIR     ${PROJECT_SOURCE_DIR}/include)

    file(GLOB headers
        ${_Header_DIR}/*.hh             ${_Header_DIR}/*.icc)

    file(GLOB sources
        ${_Source_DIR}/*.cc)

    # avoid conflicting format targets
    set(FORMAT_NAME format)
    if(TARGET format)
        set(FORMAT_NAME format-ptl)
    endif()

    add_custom_target(${FORMAT_NAME}
        COMMAND ${CLANG_FORMATTER} -i ${headers} ${sources}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Running '${CLANG_FORMATTER}' on '${_Source_DIR}'..."
        SOURCES ${headers} ${sources})

endif()
