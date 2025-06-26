#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(DICOM01)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(dicom1 dicom1.cc)
target_link_libraries(dicom1 PRIVATE G4DicomCore::G4DicomCore)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build DICOM. This is so that we can run the executable directly because it
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

# new DICOM data (in share directory)
set(DICOM_SHARE
  AltData.dat SixSlice.dat
  IM-0003-0001.dcm IM-0003-0003.dcm IM-0003-0005.dcm IM-0003-0007.dcm IM-0003-0009.dcm
  IM-0003-0002.dcm IM-0003-0004.dcm IM-0003-0006.dcm IM-0003-0008.dcm IM-0003-0010.dcm)

# copy over scripts
foreach(_script ${DICOM_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY)
endforeach()

# copy either Data.dat.old or Data.dat.new to Data.dat based on build settings
if(G4DICOM_USE_DCMTK)
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

