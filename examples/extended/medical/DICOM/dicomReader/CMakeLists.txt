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
# Install the library to 'lib' directory under CMAKE_INSTALL_PREFIX
#
set(dicomReader_DIR lib/Geant4-${Geant4_VERSION} CACHE PATH "dicomReader installation")

install(FILES ${headers} DESTINATION include/dicomReader)
install(TARGETS ${${PROJECT_NAME}_INSTALL_LIBRARIES}
    DESTINATION lib EXPORT dicomReaderTargets)
install(EXPORT dicomReaderTargets DESTINATION ${dicomReader_DIR})

set(INCLUDE_INSTALL_DIR include/dicomReader)
set(LIB_INSTALL_DIR lib)
set(CONFIG_INSTALL_DIR ${dicomReader_DIR})

include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/dicomReaderConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/dicomReaderConfig.cmake
    INSTALL_DESTINATION ${dicomReader_DIR}
    PATH_VARS INCLUDE_INSTALL_DIR CONFIG_INSTALL_DIR)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/dicomReaderConfigVersion.cmake
    VERSION ${Geant4_VERSION}
    COMPATIBILITY SameMajorVersion )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/dicomReaderConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/dicomReaderConfigVersion.cmake
    DESTINATION ${dicomReader_DIR} )
