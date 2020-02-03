#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(dicomReader)
# for shared library
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
#set(CMAKE_MACOSX_RPATH OFF)

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
include("${PROJECT_SOURCE_DIR}/cmake/DICOMUtilities.cmake")

#----------------------------------------------------------------------------
# Find DCMTK installation
#
find_package(DCMTK REQUIRED)
add_definitions(-DG4_DCMTK)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${DCMTK_INCLUDE_DIRS})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the library
#
dicom_build_library(
    BUILD_SHARED ${Geant4_shared_FOUND}
    BUILD_STATIC ${Geant4_static_FOUND}
    OUTPUT_NAME dicomReader
    TARGET_NAME dicomReader
    SOURCES ${headers} ${sources}
    LINK_LIBRARIES ${Geant4_LIBRARIES} ${DCMTK_LIBRARIES})

#----------------------------------------------------------------------------
# Configuration
#
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

#----------------------------------------------------------------------------
# Configuration for build tree
#
export(TARGETS ${${PROJECT_NAME}_INSTALL_LIBRARIES}
    FILE ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Build.cmake)

set(PROJECT_TARGETS_FILE ${PROJECT_NAME}Build.cmake)
set(PACKAGE_INCLUDE_INSTALL_DIR ${PROJECT_SOURCE_DIR}/include)
set(PACKAGE_INIT
"macro(set_and_check _var _file)
    set(\${_var} \"\${_file}\")
    if(NOT EXISTS \"\${_file}\")
        message(FATAL_ERROR \"File or directory \${_file} referenced by variable \${_var} does not exist !\")
    endif()
endmacro()

macro(check_required_components _NAME)
    foreach(comp \${\${_NAME}_FIND_COMPONENTS})
        if(NOT \${_NAME}_\${comp}_FOUND)
            if(\${_NAME}_FIND_REQUIRED_\${comp})
                set(\${_NAME}_FOUND FALSE)
            endif()
        endif()
    endforeach()
endmacro()
")

configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake @ONLY)

write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    VERSION ${Geant4_VERSION}
    COMPATIBILITY SameMajorVersion )

unset(PACKAGE_INIT)

#----------------------------------------------------------------------------
# Install the library to 'lib' directory under CMAKE_INSTALL_PREFIX
#
set(${PROJECT_NAME}_DIR ${CMAKE_INSTALL_LIBDIR}/Geant4-${Geant4_VERSION} CACHE PATH "${PROJECT_NAME} installation")

install(FILES ${headers} DESTINATION include/${PROJECT_NAME})
install(TARGETS ${${PROJECT_NAME}_INSTALL_LIBRARIES}
    DESTINATION ${CMAKE_INSTALL_LIBDIR} EXPORT ${PROJECT_NAME}Targets)
install(EXPORT ${PROJECT_NAME}Targets DESTINATION ${${PROJECT_NAME}_DIR})

set(PROJECT_TARGETS_FILE ${PROJECT_NAME}Targets.cmake)
set(INCLUDE_INSTALL_DIR include/${PROJECT_NAME})
set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})

configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_BINARY_DIR}/InstallTreeFiles/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${${PROJECT_NAME}_DIR}
    PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR)

write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/InstallTreeFiles/${PROJECT_NAME}ConfigVersion.cmake
    VERSION ${Geant4_VERSION}
    COMPATIBILITY SameMajorVersion )

install(FILES ${CMAKE_BINARY_DIR}/InstallTreeFiles/${PROJECT_NAME}Config.cmake
    ${CMAKE_BINARY_DIR}/InstallTreeFiles/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION ${${PROJECT_NAME}_DIR} )
