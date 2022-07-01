# -------------------------------------------------------------------------------------- #
# cmake settings
#
cmake_minimum_required(VERSION 3.8...3.23)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# Warn about an in-source build, but don't error out
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(STATUS "Warning! Building from the source directory is not recommended")
    message(STATUS "If unintented, please remove 'CMakeCache.txt' and 'CMakeFiles'")
    message(STATUS "and build from a separate directory")
    message(WARNING "In-source build")
endif()

# Check if project is being used directly or via add_subdirectory
set(PTL_MASTER_PROJECT ON)
if(NOT CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(PTL_MASTER_PROJECT OFF)
else()
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# Version
file(READ "${CMAKE_CURRENT_LIST_DIR}/VERSION" VERSION_STRING LIMIT_COUNT 1)
string(REGEX REPLACE "(\n|\r)" "" VERSION_STRING "${VERSION_STRING}")
string(REGEX REPLACE "[A-Za-z].*" "" VERSION_STRING "${VERSION_STRING}")
set(PTL_VERSION "${VERSION_STRING}")

# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
# project setup
#
project(
    PTL
    LANGUAGES C CXX
    VERSION ${PTL_VERSION})

# Postprocess version info to create variables for export to Version.hh
set(PTL_VERSION_STRING "${PTL_VERSION_MAJOR}.${PTL_VERSION_MINOR}.${PTL_VERSION_PATCH}")
math(EXPR PTL_VERSION_CODE
     "10000 * ${PTL_VERSION_MAJOR} + 100 * ${PTL_VERSION_MINOR} + ${PTL_VERSION_PATCH}")

# Project-local CMake settings
set(CMAKE_DIRECTORY_LABELS "PTL")
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME Development)
set(CMAKE_INSTALL_MESSAGE LAZY)
set(CMAKE_LINK_DEPENDS_NO_SHARED ON)
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_LIST_DIR}/cmake/Modules)
include(PTLCMakeUtilities)

# -------------------------------------------------------------------------------------- #
# Build/Install settings and options
include(PTLInstallDirs)
include(PTLBuildSettings)

# -------------------------------------------------------------------------------------- #
# User options
ptl_add_option(PTL_USE_TBB "Enable TBB" ON)
ptl_add_option(PTL_USE_LOCKS "Enable mutex locking in task subqueues for extra safety"
               OFF)
ptl_add_option(PTL_INSTALL_HEADERS "Install the headers" ON)
ptl_add_option(PTL_INSTALL_CONFIG "Install the cmake configuration" ON)

if(DEFINED PTL_DEVELOPER_INSTALL)
    set(PTL_INSTALL_HEADERS
        ${PTL_DEVELOPER_INSTALL}
        CACHE BOOL "Set via PTL_DEVELOPER_INSTALL" FORCE)
    set(PTL_INSTALL_CONFIG
        ${PTL_DEVELOPER_INSTALL}
        CACHE BOOL "Set via PTL_DEVELOPER_INSTALL" FORCE)
endif()

# -------------------------------------------------------------------------------------- #
# Build Dependencies - Threads
if(NOT WIN32)
    set(CMAKE_THREAD_PREFER_PTHREAD ON)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
endif()
find_package(Threads REQUIRED)

# * TBB
if(PTL_USE_TBB)
    find_package(TBB 2017 REQUIRED)
endif()

# -------------------------------------------------------------------------------------- #
# PTL Primary Build
add_subdirectory(src)

# -------------------------------------------------------------------------------------- #
# CMake/PkgConfig Support files
include(PTLPackageConfigHelpers)

# -------------------------------------------------------------------------------------- #
# Examples build/test
if(EXISTS ${CMAKE_CURRENT_LIST_DIR}/examples)
    ptl_add_option(PTL_BUILD_EXAMPLES "Build examples" OFF)
    if(PTL_BUILD_EXAMPLES)
        set(PTL_DIR ${CMAKE_BINARY_DIR})
        add_subdirectory(examples)
    endif()
endif()

# -------------------------------------------------------------------------------------- #
# Reporting if master project
if(PTL_MASTER_PROJECT)
    ptl_print_features()
endif()
