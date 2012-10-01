cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test69)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test69 test69.cc src/*.cc)

file(GLOB _MACRO_FILES ${CMAKE_CURRENT_SOURCE_DIR}/macros/*.in)

#---Test definitions----------------------------------------------------------------------------------
# We define a test69-build test that only builds the executable.
# We also define a normal test for each existing .in file.
GEANT4_ADD_TEST(test69-build BUILD test69 LABELS Nightly Continuous) 

foreach(_macro_file ${_MACRO_FILES})
  get_filename_component(_test_extension ${_macro_file} NAME_WE)

  GEANT4_ADD_TEST(test69-${_test_extension}
                   COMMAND test69 ${_macro_file}
                   DEPENDS test69-build
                   LABELS Nightly Continuous
                   ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

endforeach(_macro_file _MACRO_FILES)

