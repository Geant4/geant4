################################################################################
#   cmake settings
#
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

# Add allocation export symbol for the PLT module
add_definitions(-DPTL_ALLOC_EXPORT)

cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0042 NEW)

# Check if project is being used directly or via add_subdirectory
set(PTL_MASTER_PROJECT ON)
if(NOT CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(PTL_MASTER_PROJECT OFF)
  set(CMAKE_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/Geant4")
endif()

################################################################################
#   version
#
file(READ "${CMAKE_CURRENT_LIST_DIR}/VERSION" VERSION_STRING LIMIT_COUNT 1)
string(REGEX REPLACE "(\n|\r)" "" VERSION_STRING "${VERSION_STRING}")
string(REGEX REPLACE "[A-Za-z].*" "" VERSION_STRING "${VERSION_STRING}")
set(PTL_VERSION "${VERSION_STRING}")

################################################################################
#   project
#
project(PTL LANGUAGES C CXX VERSION ${PTL_VERSION})
if(NOT PTL_MASTER_PROJECT)
    unset(${PROJECT_NAME}_C_FLAGS CACHE)
    unset(${PROJECT_NAME}_CXX_FLAGS CACHE)
endif()
set(CMAKE_DIRECTORY_LABELS "PTL")

################################################################################
#
#
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})
set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME Development)

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "CMake build type" FORCE)
endif()

################################################################################

include(MacroUtilities)
#include(Options)
include(Compilers)
include(GNUInstallDirs)
include(ProjectSettings)
include(BuildSettings)
include(Packages)
include(ClangFormat)

################################################################################
#   PTL source
################################################################################

cmake_policy(SET CMP0053 NEW)

################################################################################
#                                                                              #
#                               General                                        #
#                                                                              #
################################################################################

if(PTL_USE_TBB)
    list(APPEND ${PROJECT_NAME}_DEFINITIONS PTL_USE_TBB)
endif()

# libraries to install
set(INSTALL_LIBRARIES )

#------------------------------------------------------------------------------#
# Locate sources and headers for this project
# - headers are included so they will show up in IDEs
file(GLOB_RECURSE ptl_headers ${CMAKE_CURRENT_LIST_DIR}/include/PTL/*.hh
                              ${CMAKE_CURRENT_LIST_DIR}/include/PTL/*.icc)
file(GLOB_RECURSE ptl_sources ${CMAKE_CURRENT_LIST_DIR}/src/*.cc)

get_property(PTL_INTERFACE_LIBRARIES GLOBAL PROPERTY ${PROJECT_NAME}_INTERFACE_LIBRARIES)
export(TARGETS ${PTL_INTERFACE_LIBRARIES} FILE ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-interface.cmake)
################################################################################
#                                                                              #
#                              PTL Library                                     #
#                                                                              #
################################################################################

if(BUILD_SHARED_LIBS)
    build_library(VERSION
        TYPE SHARED
        TARGET_NAME ptl-shared
        OUTPUT_NAME G4ptl
        SOURCES ${ptl_headers} ${ptl_sources}
        EXTRA_ARGS POSITION_INDEPENDENT_CODE ON
                   WINDOWS_EXPORT_ALL_SYMBOLS ON)

    target_link_libraries(ptl-shared PUBLIC ptl-external-packages)
    target_link_libraries(ptl-shared PRIVATE ptl-external-libraries ptl-compile-options)

    target_include_directories(ptl-shared PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

    list(APPEND INSTALL_LIBRARIES ptl-shared)
    export(TARGETS ptl-shared FILE ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-shared.cmake)
endif()

if(BUILD_STATIC_LIBS)

    build_library(
        TYPE STATIC
        TARGET_NAME ptl-static
        OUTPUT_NAME $<IF:$<PLATFORM_ID:Windows>,G4ptl-static,G4ptl>
        SOURCES ${ptl_headers} ${ptl_sources}
        EXTRA_ARGS POSITION_INDEPENDENT_CODE ON)

    target_link_libraries(ptl-static PUBLIC ptl-external-packages)
    target_link_libraries(ptl-static PRIVATE ptl-external-libraries ptl-compile-options)

    target_include_directories(ptl-static PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

    list(APPEND INSTALL_LIBRARIES ptl-static)
    export(TARGETS ptl-static FILE ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-static.cmake)
endif()


################################################################################
#                                                                              #
#                              Installation                                    #
#                                                                              #
################################################################################

# Install the targets and export libraries
install(TARGETS ${INSTALL_LIBRARIES}
    EXPORT ${PROJECT_NAME}Targets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT Development
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT Runtime
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT Runtime)

# install export
install(EXPORT ${PROJECT_NAME}Targets
    NAMESPACE PTL::
    DESTINATION ${CMAKE_INSTALL_CONFIGDIR}
    COMPONENT Development)

# headers
install(FILES ${ptl_headers}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/PTL
    COMPONENT Development)

# set
set(PTL_BUILT_LIBRARIES "${INSTALL_LIBRARIES}" CACHE INTERNAL "Built libraries" FORCE)

################################################################################
#   Installation and info
################################################################################

include(ConfigurePackage)

if(PTL_MASTER_PROJECT)
    print_features()
endif()

