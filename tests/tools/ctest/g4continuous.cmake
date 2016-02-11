cmake_minimum_required(VERSION 2.6)

#---Common Geant4 CTest script----------------------------------------------
include(${CTEST_SCRIPT_DIRECTORY}/g4common.cmake)

#---mark continuous in build name-------------------------------------------
set(CTEST_BUILD_NAME c_${CTEST_BUILD_NAME})

#---Addional CTest settings-------------------------------------------------
set(CTEST_UPDATE_OPTIONS "${CTEST_UPDATE_OPTIONS} -p")  # Add proposed tags

#---Clean or not clean the Binary-------------------------------------------
GET_DATE(date)
if(NOT EXISTS ${CTEST_BINARY_DIRECTORY}/${date} )
  ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
endif()

#---CTest commands----------------------------------------------------------
ctest_start (Continuous)
ctest_update()
ctest_configure(BUILD   ${CTEST_BINARY_DIRECTORY}
                SOURCE  ${CTEST_SOURCE_DIRECTORY}
                OPTIONS "${CTEST_CONFIG_OPTIONS}")
ctest_build(BUILD ${CTEST_BINARY_DIRECTORY})
ctest_test(PARALLEL_LEVEL ${ncpu}
               INCLUDE_LABEL "Nightly|Continuous"
               EXCLUDE "largeN$"
               INCLUDE "^test|^validate")
ctest_submit()

file(WRITE ${CTEST_BINARY_DIRECTORY}/${date} "timestamp file")
