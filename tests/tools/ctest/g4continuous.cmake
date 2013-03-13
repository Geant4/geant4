cmake_minimum_required(VERSION 2.6)

#---Common Geant4 CTest script----------------------------------------------
include(${CTEST_SCRIPT_DIRECTORY}/g4common.cmake)

#---Addional CTest settings-------------------------------------------------
set(CTEST_UPDATE_OPTIONS "${CTEST_UPDATE_OPTIONS} -p")  # Add proposed tags

#---CTest commands----------------------------------------------------------
ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
GET_TIME(current_time)
while(${current_time} LESS 2300)
  ctest_start (Continuous)
  ctest_update(RETURN_VALUE updates)
  if(updates GREATER 0)
    ctest_configure(BUILD   ${CTEST_BINARY_DIRECTORY} 
                    SOURCE  ${CTEST_SOURCE_DIRECTORY}
                    OPTIONS "${CTEST_CONFIG_OPTIONS}")
    ctest_build(BUILD ${CTEST_BINARY_DIRECTORY})
    ctest_test(PARALLEL_LEVEL ${ncpu}
               INCLUDE_LABEL "Nightly|Continuous"
               EXCLUDE "largeN$"
               INCLUDE "^test|^validate")
    ctest_submit()
  endif()
  ctest_sleep(600)
  GET_TIME(current_time)  
endwhile()
