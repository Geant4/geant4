
#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(DICOM2)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# append CMAKE_PREFIX_PATH to help find DICOM in recommended installation
# location (${Geant4_DIR}) and the location found when testing
if(NOT DEFINED DICOM_DIR)
    # recommended installation place
    list(APPEND CMAKE_PREFIX_PATH ${Geant4_DIR})
    # found here when internally testing
    list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_BINARY_DIR}/../DICOM")
endif()

# find DICOM package
find_package(DICOM REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include ${DICOM_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(dicom2 dicom2.cc ${sources} ${headers})
target_link_libraries(dicom2 ${Geant4_LIBRARIES} ${DICOM_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
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

foreach(_script ${DICOM_SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/data/${_script}
        ${PROJECT_BINARY_DIR}/${_script}
        COPYONLY)
endforeach()

if(DICOM_USE_DCMTK)
    configure_file(${PROJECT_SOURCE_DIR}/data/Data.dat.new
        ${PROJECT_BINARY_DIR}/Data.dat COPYONLY)
else()
    configure_file(${PROJECT_SOURCE_DIR}/data/Data.dat.old
        ${PROJECT_BINARY_DIR}/Data.dat COPYONLY)
endif()

foreach(_script ${DICOM_MACROS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${PROJECT_BINARY_DIR}/${_script}
        COPYONLY)
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS dicom2 DESTINATION bin)
