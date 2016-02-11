cmake_minimum_required(VERSION 2.8)

#---Common Geant4 CTest script----------------------------------------------
include(${CTEST_SCRIPT_DIRECTORY}/g4common.cmake)

#---Addional CTest settings-------------------------------------------------
#set(CTEST_UPDATE_OPTIONS "${CTEST_UPDATE_OPTIONS} -p")  # Add proposed tags

#---mark experimental in build name-----------------------------------------
set(CTEST_BUILD_NAME e_${CTEST_BUILD_NAME})

#---CTest commands----------------------------------------------------------
ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
ctest_start("Experimental")
ctest_update()
ctest_configure(BUILD   ${CTEST_BINARY_DIRECTORY} 
                SOURCE  ${CTEST_SOURCE_DIRECTORY}
                OPTIONS "${CTEST_CONFIG_OPTIONS}")
ctest_build(BUILD ${CTEST_BINARY_DIRECTORY})
#ctest_test(PARALLEL_LEVEL ${ncpu} INCLUDE_LABEL "Experimental")
ctest_test(PARALLEL_LEVEL ${ncpu} INCLUDE_LABEL "Nightly")
ctest_submit()




