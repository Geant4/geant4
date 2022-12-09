# PTL Build Settings and Options
include(PTLCMakeUtilities)
include(CheckCXXCompilerFlag)

# -------------------------------------------------------------------------------------- #
#
# Build Type(s)
#
# -------------------------------------------------------------------------------------- #
set(__ptl_default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(
        STATUS
            "Forcing PTL CMAKE_BUILD_TYPE to '${__ptl_default_build_type}' as none was specified"
        )
    set(CMAKE_BUILD_TYPE
        "${__ptl_default_build_type}"
        CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel"
                                                 "RelWithDebInfo")
endif()

# -------------------------------------------------------------------------------------- #
#
# Output Directories for build products
#
# -------------------------------------------------------------------------------------- #
# If master project, set the output directory (critical on Windows and Xcode) Otherwise we
# leave this at the discretion of the parent project
#
if(PTL_MASTER_PROJECT)
    set(_BIN_DIR ${PROJECT_BINARY_DIR})
    if(WIN32)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${_BIN_DIR}/outputs/runtime)
        set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${_BIN_DIR}/outputs/library)
        set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${_BIN_DIR}/outputs/archive)
    else()
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${_BIN_DIR})
        set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${_BIN_DIR})
        set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${_BIN_DIR})
    endif()
endif()

# -------------------------------------------------------------------------------------- #
#
# Library types
#
# -------------------------------------------------------------------------------------- #
ptl_add_option(BUILD_STATIC_LIBS "Build static library" ON)
ptl_add_option(BUILD_SHARED_LIBS "Build shared library" ON)

if(NOT "${CMAKE_PROJECT_NAME}" STREQUAL "${PROJECT_NAME}")
    ptl_add_option(BUILD_OBJECT_LIBS "Build object library (only valid as subproject)"
                   OFF)
    set(PTL_EXCLUDE_FROM_ALL EXCLUDE_FROM_ALL)
    set(PTL_BUILD_LIBS_ERROR_DESC
        "BUILD_STATIC_LIBS, BUILD_SHARED_LIBS, and BUILD_OBJECT_LIBS are all OFF")
else()
    set(BUILD_OBJECT_LIBS OFF)
    set(PTL_EXCLUDE_FROM_ALL)
    set(PTL_BUILD_LIBS_ERROR_DESC
        "Neither BUILD_STATIC_LIBS nor BUILD_SHARED_LIBS are set")
endif()

if((NOT BUILD_SHARED_LIBS)
   AND (NOT BUILD_STATIC_LIBS)
   AND (NOT BUILD_OBJECT_LIBS))
    message(FATAL_ERROR "${PTL_BUILD_LIBS_ERROR_DESC}. One must be ON")
endif()

# -------------------------------------------------------------------------------------- #
#
# C/C++ Standard and Flags
#
# -------------------------------------------------------------------------------------- #
if(WIN32)
    set(CMAKE_CXX_STANDARD
        14
        CACHE STRING "C++ STL standard")
else()
    set(CMAKE_CXX_STANDARD
        11
        CACHE STRING "C++ STL standard")
endif()
ptl_message_on_change(CMAKE_CXX_STANDARD "Building PTL with CMAKE_CXX_STANDARD")

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add allocation export symbol for the PTL module Should this be a
# target_compile_definition?
add_compile_definitions(PTL_ALLOC_EXPORT)

# -------------------------------------------------------------------------------------- #
#
# Testing/Profiling Options
#
# -------------------------------------------------------------------------------------- #
# Coverage compile/link settings See also:
# https://discourse.cmake.org/t/best-practice-for-multi-config-build-options/2049
ptl_add_option(PTL_USE_COVERAGE "Enable code coverage" OFF)
if(PTL_USE_COVERAGE)
    if((NOT CMAKE_BUILD_TYPE STREQUAL "Debug") AND (NOT CMAKE_CONFIGURATION_TYPES))
        message(
            WARNING
                "PTL_USE_COVERAGE is ON, but the build mode is not 'Debug': coverage testing will be disabled."
            )
        set(PTL_USE_COVERAGE
            OFF
            CACHE BOOL "Enable code coverage (disabled)" FORCE)
    else()
        include(PTLCoverage)
    endif()
endif()

# Sanitizer compile/link flags See also:
# https://discourse.cmake.org/t/best-practice-for-multi-config-build-options/2049
ptl_add_option(PTL_USE_SANITIZER "Enable -fsanitize=<type>" OFF)
if(PTL_USE_SANITIZER)
    set(__ptl_sanitizer_type_default "thread")
    if(PTL_SANITIZER_TYPE)
        set(__ptl_sanitizer_type_default "${PTL_SANITIZER_TYPE}")
    endif()

    set(PTL_SANITIZER_TYPE
        "${__ptl_sanitizer_type_default}"
        CACHE STRING "Sanitizer type (-fsanitize=<type>)")
    set_property(CACHE PTL_SANITIZER_TYPE PROPERTY STRINGS "thread" "address" "undefined")
    ptl_message_on_change(PTL_SANITIZER_TYPE "Building PTL with sanitizer type")

    include(PTLSanitize)
endif()

# Clang Tooling
if(PTL_MASTER_PROJECT)
    include(PTLFormatting)
endif()

ptl_add_option(PTL_USE_CLANG_TIDY "Enable running clang-tidy on sources" OFF)
if(PTL_USE_CLANG_TIDY)
    include(PTLLinting)
endif()
