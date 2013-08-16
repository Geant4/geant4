#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(DICOM)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
                    
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(DICOM DICOM.cc ${headers} ${sources})

target_link_libraries(DICOM ${Geant4_LIBRARIES} )


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build DICOM. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

# the macros
set(DICOM_MACROS
    run.mac vis.mac
    )
 
# original set of DICOM data
set(DICOM_SCRIPTS
    1.dcm 2.dcm 3.dcm 
    1.g4  2.g4 3.g4 
    1.g4dcm 2.g4dcm 3.g4dcm
    ColourMap.dat CT2Density.dat 
    Data.dat 
  )
  
# new DICOM data (in share directory)
set(DICOM_SHARE
    AltData.dat SixSlice.dat
    IM-0003-0001.dcm IM-0003-0003.dcm IM-0003-0005.dcm IM-0003-0007.dcm IM-0003-0009.dcm
    IM-0003-0002.dcm IM-0003-0004.dcm IM-0003-0006.dcm IM-0003-0008.dcm IM-0003-0010.dcm
  )

foreach(_script ${DICOM_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

foreach(_script ${DICOM_SHARE})
  configure_file(
    ${PROJECT_SOURCE_DIR}/share/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

foreach(_script ${DICOM_MACROS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS DICOM DESTINATION bin)

