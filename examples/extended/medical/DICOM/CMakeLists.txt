#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(G4DICOM)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# See the documentation for a guide on how to enable/disable specific components
#
find_package(Geant4 REQUIRED ui_all vis_all)

#----------------------------------------------------------------------------
# We can optional build these examples with support for the DCMTK libraries,
# and use of HEAD data.
option(G4DICOM_USE_DCMTK "Build with DCMTK support" OFF)
option(G4DICOM_USE_HEAD "Download and use HEAD data" OFF)

# For HEAD data, we need to download the data
#----------------------------------------------------------------------------
# Download DICOM_HEAD data
#
# enable option if environment variable set (backwards-compat)
if(G4DICOM_USE_HEAD)
  message(STATUS "Enabling DICOM_HEAD data download")
  set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
  include(DownloadDICOMData)
endif()

# If we choose to use DCMTK, we must find it, then build the interface to it 
# TODO: Make this work with DCMTK 3
if(G4DICOM_USE_DCMTK)
  add_subdirectory(G4DicomReader)
endif()

# Build library of functionality shared between the two example applications
add_subdirectory(G4DicomCore)

# Build the example applications
add_subdirectory(DICOM1)
add_subdirectory(DICOM2)


