
#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(DICOM2)

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
target_include_directories(dicom2 PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(dicom2 PRIVATE G4DicomCore::G4DicomCore ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
# the macros
set(DICOM_MACROS run.mac vis.mac)

# original set of DICOM data
set(DICOM_SCRIPTS
  1.dcm   2.dcm   3.dcm
  1.g4    2.g4    3.g4
  1.g4dcm 2.g4dcm 3.g4dcm
  ColourMap.dat CT2Density.dat
  Data.dat.new Data.dat.old
  Data.dat.new_dens Data.partial.dat)

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
