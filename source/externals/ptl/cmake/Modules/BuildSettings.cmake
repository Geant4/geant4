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
ptl_add_interface_library(ptl-external-libraries)

# ---------------------------------------------------------------------------- #
#
#set(SANITIZE_TYPE leak CACHE STRING "-fsantitize=<TYPE>")

# ---------------------------------------------------------------------------- #
# if master project, set the output directory (critical on Windows and Xcode)
#
if("${CMAKE_PROJECT_NAME}" STREQUAL "${PROJECT_NAME}")
    set(_BIN_DIR ${PROJECT_BINARY_DIR})
    if(WIN32)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${_BIN_DIR}/outputs/Runtime)
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

