cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(volumes_tests)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})


file(GLOB tests_files ${PROJECT_SOURCE_DIR}/*.cc)


foreach(_tests ${tests_files})
 GET_FILENAME_COMPONENT(tests_files_we ${_tests} NAME_WE)
 add_executable(${tests_files_we} ${_tests})
 target_link_libraries(${tests_files_we} ${Geant4_LIBRARIES})
endforeach()


