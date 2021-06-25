################################################################################
#
#        Handles the build settings
#
################################################################################
#
set(LIBNAME ptl)

include(GNUInstallDirs)
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
include(Compilers)
include(MacroUtilities)

ptl_add_interface_library(ptl-compile-options)
ptl_add_interface_library(ptl-public-options)
ptl_add_interface_library(ptl-external-libraries)
ptl_add_interface_library(ptl-sanitizer-options)

# ---------------------------------------------------------------------------- #
#
#set(SANITIZE_TYPE leak CACHE STRING "-fsantitize=<TYPE>")

# ---------------------------------------------------------------------------- #
# if master project, set the output directory (critical on Windows and Xcode)
#
if("${CMAKE_PROJECT_NAME}" STREQUAL "${PROJECT_NAME}")
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

# ---------------------------------------------------------------------------- #
#  debug macro
#
string(TOUPPER "${CMAKE_BUILD_TYPE}" UC_BUILD_TYPE)
if("${UC_BUILD_TYPE}" STREQUAL "DEBUG")
    target_compile_definitions(ptl-compile-options INTERFACE DEBUG)
else()
    target_compile_definitions(ptl-compile-options INTERFACE NDEBUG)
endif()

if(PTL_USE_LOCKS)
    target_compile_definitions(ptl-public-options INTERFACE PTL_USE_LOCKS)
endif()

if(PTL_USE_COVERAGE)
    target_compile_options(ptl-public-options INTERFACE
        $<BUILD_INTERFACE:-fprofile-arcs>
        $<BUILD_INTERFACE:-ftest-coverage>)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(ptl-public-options INTERFACE
            $<BUILD_INTERFACE:-fprofile-abs-path>
            $<BUILD_INTERFACE:--coverage>)
        set_target_properties(ptl-public-options PROPERTIES
            INTERFACE_LINK_OPTIONS
                $<BUILD_INTERFACE:--coverage>)
    else()
        set_target_properties(ptl-public-options PROPERTIES
            INTERFACE_LINK_OPTIONS
                $<BUILD_INTERFACE:-fprofile-arcs>)
    endif()
endif()

if(PTL_USE_SANITIZER AND PTL_SANITIZER_TYPE)
    if("${PTL_SANITIZER_TYPE}" STREQUAL "thread")
        target_compile_options(ptl-sanitizer-options INTERFACE
            -fsanitize=${PTL_SANITIZER_TYPE})
    else()
        target_compile_options(ptl-sanitizer-options INTERFACE
            -fsanitize=${PTL_SANITIZER_TYPE}
            -fno-optimize-sibling-calls
            -fno-omit-frame-pointer
            -fno-inline-functions)
    endif()
    target_link_options(ptl-sanitizer-options INTERFACE
        -fsanitize=${PTL_SANITIZER_TYPE})
endif()

