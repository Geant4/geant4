# * Clang-format only for master project
if(PTL_CLANG_FORMATTER MATCHES ".*-6")
    unset(PTL_CLANG_FORMATTER CACHE)
endif()

find_program(PTL_CLANG_FORMATTER NAMES clang-format-9 clang-format-mp-9.0 clang-format)
mark_as_advanced(PTL_CLANG_FORMATTER)

find_program(PTL_CMAKE_FORMATTER NAMES cmake-format)
mark_as_advanced(PTL_CMAKE_FORMATTER)

set(PTL_FORMAT_TARGET format)
if(TARGET format)
    set(PTL_FORMAT_TARGET format-ptl)
endif()

if(PTL_CLANG_FORMATTER)
    set(_Source_DIR ${PROJECT_SOURCE_DIR}/source)
    set(_Example_DIR ${PROJECT_SOURCE_DIR}/examples)

    file(GLOB_RECURSE headers ${_Source_DIR}/*.hh ${_Source_DIR}/*.icc
         ${_Example_DIR}/*.hh ${_Example_DIR}/*.h)

    file(GLOB_RECURSE sources ${_Source_DIR}/*.cc ${_Source_DIR}/*.c ${_Example_DIR}/*.cc
         ${_Example_DIR}/*.cu)

    add_custom_target(
        ${PTL_FORMAT_TARGET}-source
        COMMAND ${PTL_CLANG_FORMATTER} -i ${headers} ${sources}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT
            "Running '${PTL_CLANG_FORMATTER}' on '${_Source_DIR}' and '${_Example_DIR}..."
        SOURCES ${headers} ${sources})
endif()

if(PTL_CMAKE_FORMATTER)
    set(_Source_DIR ${PROJECT_SOURCE_DIR}/source)
    set(_Example_DIR ${PROJECT_SOURCE_DIR}/examples)

    file(GLOB_RECURSE cmake_files "${PROJECT_SOURCE_DIR}/**/CMakeLists.txt"
         "${PROJECT_SOURCE_DIR}/*.cmake" "${PROJECT_SOURCE_DIR}/*.cmake.in")
    list(REMOVE_ITEM cmake_files "${PROJECT_SOURCE_DIR}/cmake/Modules/FindTBB.cmake")

    add_custom_target(
        ${PTL_FORMAT_TARGET}-cmake
        COMMAND ${PTL_CMAKE_FORMATTER} -i ${cmake_files}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Running '${PTL_CMAKE_FORMATTER}'..."
        SOURCES ${cmake_files})
endif()

foreach(_FORMAT_TARGET ${PTL_FORMAT_TARGET}-source)
    # avoid conflicting format targets
    if(NOT TARGET ${PTL_FORMAT_TARGET})
        add_custom_target(${PTL_FORMAT_TARGET})
    endif()
    if(TARGET ${_FORMAT_TARGET})
        add_dependencies(${PTL_FORMAT_TARGET} ${_FORMAT_TARGET})
    endif()
endforeach()
