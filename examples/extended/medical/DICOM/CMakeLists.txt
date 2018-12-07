#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(DICOM)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# DICOM configure options
#
include(CMakeDependentOption)
# this macro checks for environment variable ${VAR} and sets ENV_${VAR} if
# if is a non-empty string. If empty string and ${VAR} is defined in CMake
# cache, then set ENV_${VAR} to it's value
# NOTE: environment variable, if set, overrides cache value
macro(check_environment VAR)
    set(ENV_${VAR} "$ENV{${VAR}}")
    if("${ENV_${VAR}}" STREQUAL "" AND DEFINED ${VAR})
        set(ENV_${VAR} ${${VAR}})
    else()
        set(ENV_${VAR} OFF) # default to off
    endif()
endmacro(check_environment VAR)

check_environment(DICOM_USE_DCMTK)
check_environment(DICOM_USE_HEAD)
# enable option if environment variable set (backwards-compat)
option(DICOM_USE_DCMTK "DICOM with DCMTK support" ${ENV_DICOM_USE_DCMTK})
option(DICOM_USE_HEAD "Download DICOM_HEAD data" ${ENV_DICOM_USE_HEAD})

if(DICOM_USE_HEAD)
  add_definitions(-DDICOM_USE_HEAD)
endif()


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
include("${PROJECT_SOURCE_DIR}/dicomReader/cmake/DICOMUtilities.cmake")

#----------------------------------------------------------------------------
# Add dicomReader subdirectory
#
if(DICOM_USE_DCMTK)
    message(STATUS "${PROJECT_NAME}: Using DCMTK")
    find_package(DCMTK REQUIRED)
    add_definitions(-DG4_DCMTK)
    add_subdirectory(dicomReader)
    set(DICOM_READER_LIBRARY dicomReader${_geant4_lib_use_suffix})
endif()

#----------------------------------------------------------------------------
# Download DICOM_HEAD data
#
# enable option if environment variable set (backwards-compat)
if(DICOM_USE_HEAD)
    message(STATUS "${PROJECT_NAME}: Enabling DICOM_HEAD data download")
    set(CMAKE_MODULE_PATH
        ${CMAKE_SOURCE_DIR}/cmake
        ${CMAKE_MODULE_PATH})
  include(DownloadDICOMData)
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${PROJECT_SOURCE_DIR}/dicomReader/include
                    ${Geant4_INCLUDE_DIR}
                    ${DCMTK_INCLUDE_DIRS})
                    
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# List any source specific properties here



#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
dicom_build_library(
    BUILD_SHARED ${Geant4_shared_FOUND}
    BUILD_STATIC ${Geant4_static_FOUND}
    OUTPUT_NAME DICOM
    TARGET_NAME DICOM-library
    SOURCES ${headers} ${sources}
    LINK_LIBRARIES ${Geant4_LIBRARIES} ${DICOM_READER_LIBRARY} ${DCMTK_LIBRARIES})

# DICOM-library-target is set in dicom_build_library
add_executable(DICOM DICOM.cc)
target_link_libraries(DICOM DICOM::library)


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build DICOM. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

# the macros
set(DICOM_MACROS run.mac vis.mac)

# original set of DICOM data
set(DICOM_SCRIPTS
    1.dcm 2.dcm 3.dcm
    1.g4  2.g4 3.g4
    1.g4dcm 2.g4dcm 3.g4dcm
    ColourMap.dat CT2Density.dat
    Data.dat.new Data.dat.old
    Data.dat.new_dens Data.partial.dat
)
 
# new DICOM data (in share directory)
set(DICOM_SHARE
    AltData.dat SixSlice.dat
    IM-0003-0001.dcm IM-0003-0003.dcm IM-0003-0005.dcm IM-0003-0007.dcm IM-0003-0009.dcm
    IM-0003-0002.dcm IM-0003-0004.dcm IM-0003-0006.dcm IM-0003-0008.dcm IM-0003-0010.dcm
)

# copy over scripts
foreach(_script ${DICOM_SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${PROJECT_BINARY_DIR}/${_script}
        COPYONLY)
endforeach()

# copy either Data.dat.old or Data.dat.new to Data.dat based on build settings
if(DICOM_USE_DCMTK)
    configure_file(${PROJECT_SOURCE_DIR}/Data.dat.new
        ${PROJECT_BINARY_DIR}/Data.dat COPYONLY)
else()
    configure_file(${PROJECT_SOURCE_DIR}/Data.dat.old
        ${PROJECT_BINARY_DIR}/Data.dat COPYONLY)
endif()

# copy over files in share
foreach(_script ${DICOM_SHARE})
    configure_file(
        ${PROJECT_SOURCE_DIR}/share/${_script}
        ${PROJECT_BINARY_DIR}/${_script}
        COPYONLY)
endforeach()

# copy over macros
foreach(_script ${DICOM_MACROS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${PROJECT_BINARY_DIR}/${_script}
        COPYONLY)
endforeach()

# ensure files have correct compile definitions
if(DICOM_USE_DCMTK)
    set_source_files_properties( ${sources}
        PROPERTIES COMPILE_DEFINITIONS G4_DCMTK)
endif()

#----------------------------------------------------------------------------
# Configuration for export and installation
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
# Install the project under CMAKE_INSTALL_PREFIX
#
# by default install the package configuration to the Geant4 installation tree
set(${PROJECT_NAME}_DIR ${CMAKE_INSTALL_LIBDIR}/Geant4-${Geant4_VERSION}
    CACHE PATH "${PROJECT_NAME} installation")

install(FILES ${headers} DESTINATION include/${PROJECT_NAME})
install(TARGETS DICOM DESTINATION bin)
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

set(MSG "DICOM settings:")
set(MSG "${MSG}\n   - DICOM_USE_DCMTK: ${DICOM_USE_DCMTK}")
set(MSG "${MSG}\n   - DICOM_USE_HEAD: ${DICOM_USE_HEAD}")
message(STATUS "${MSG}")
