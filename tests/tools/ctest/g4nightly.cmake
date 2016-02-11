cmake_minimum_required(VERSION 2.8)

#---Common Geant4 CTest script----------------------------------------------
include(${CTEST_SCRIPT_DIRECTORY}/g4common.cmake)

#---CTest commands----------------------------------------------------------
ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
ctest_start("Nightly")
ctest_update()
ctest_configure(BUILD   ${CTEST_BINARY_DIRECTORY} 
                SOURCE  ${CTEST_SOURCE_DIRECTORY}
                OPTIONS "${CTEST_CONFIG_OPTIONS}")
ctest_build(BUILD ${CTEST_BINARY_DIRECTORY})
ctest_test(PARALLEL_LEVEL ${ncpu} INCLUDE_LABEL "Nightly")
ctest_submit()
