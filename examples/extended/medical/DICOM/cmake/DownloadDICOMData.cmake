# Script to download DICOM_HEAD data for DICOM example
# called from the DICOM example CMakeLists.txt
#
# The script will download the G4DICOM.1.0.tar.gz from the Geant4 Download site
# and unpack it in the example build directory.
# An error is issued if the download or unpack operation fails.
#
# Inspired by test65/DownloadLENDData.cmake

set(DICOMDATA_FILENAME "G4DICOM.1.1.tar.gz")
set(DICOMDATA_WORKING_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(DICOMDATA_LOCAL_FILENAME "${DICOMDATA_WORKING_DIR}/${DICOMDATA_FILENAME}")
set(DICOMDATA_LOCAL_ROOTDIR "${DICOMDATA_WORKING_DIR}/DICOM1.1")
set(DICOMDATA_INSTALL_DIR "${Geant4_DATASET_G4ENSDFSTATE_PATH}")
get_filename_component(DICOMDATA_INSTALL_DIR "${DICOMDATA_INSTALL_DIR}" PATH)

# The Geant4 Download site
set(DICOM_DATA_URL "http://cern.ch/geant4/support/source/${DICOMDATA_FILENAME}")
message(STATUS "DICOM_DATA_URL ${DICOM_DATA_URL}")

# Return if data directory is already present
set(DICOMDATA_NEEDS_DOWNLOAD TRUE)
if (EXISTS "${DICOMDATA_LOCAL_ROOTDIR}")
  message(STATUS "DICOM example: DICOM_HEAD data found, skipping download")
  message(STATUS "Installing '${DICOMDATA_LOCAL_ROOTDIR}' to '${DICOMDATA_INSTALL_DIR}'...")

  install(DIRECTORY ${DICOMDATA_LOCAL_ROOTDIR}
      DESTINATION ${DICOMDATA_INSTALL_DIR})
  return()
endif()

# Download tar file
message(STATUS "DICOM example: DICOM_HEAD data not found, going to download")
file(DOWNLOAD "${DICOM_DATA_URL}" "${DICOMDATA_LOCAL_FILENAME}"
  SHOW_PROGRESS
  INACTIVITY_TIMEOUT 1200
  TIMEOUT 3000
  STATUS DownloadStatus
)
if(DownloadReturnStatus)
  message(ERROR  "DICOM example: download data FAILED: ${DownloadReturnStatus}, ${DownloadStringReturnStatus}")
  return()
else()
  message(STATUS "DICOM example: download data OK")
endif()

# Unpack downloaded tarball
message(STATUS "Going to unpack: ${DICOMDATA_LOCAL_FILENAME}")
execute_process(
  COMMAND ${CMAKE_COMMAND} -E tar xfz "${DICOMDATA_LOCAL_FILENAME}"
  WORKING_DIRECTORY ${DICOMDATA_WORKING_DIR}
  OUTPUT_QUIET
  RESULT_VARIABLE __dicomdata_untar_result
)

if(__dicomdata_untar_result)
  message(ERROR  "DICOM example: failed to untar file: ${DICOMDATA_LOCAL_FILENAME}")
else()
  message(STATUS "DICOM example: untarred '${DICOMDATA_LOCAL_FILENAME}' OK")
endif()

message(STATUS "Installing '${DICOMDATA_LOCAL_ROOTDIR}' to '${DICOMDATA_INSTALL_DIR}'...")

install(DIRECTORY ${DICOMDATA_LOCAL_ROOTDIR}
    DESTINATION ${DICOMDATA_INSTALL_DIR})
