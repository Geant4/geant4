# - Basic setup for testing Geant4 using CMake/CTest
#
#-----------------------------------------------------------------------
# Integration and unit tests
# - "ENABLE_TESTING" means all tests under tests/
option(GEANT4_ENABLE_TESTING "Enable and define all the tests of the project" OFF)
geant4_add_feature(GEANT4_ENABLE_TESTING "Enable and define all the tests of the project")
mark_as_advanced(GEANT4_ENABLE_TESTING)

# - "BUILD_TESTS" means all 'tests' in individual categories.
option(GEANT4_BUILD_TESTS "Build all the tests of the project" OFF)
geant4_add_feature(GEANT4_BUILD_TESTS "Build all the tests of the project")
mark_as_advanced(GEANT4_BUILD_TESTS)

#-----------------------------------------------------------------------
# Configure CTest and relevant Geant4 settings, if required
#
if(GEANT4_ENABLE_TESTING)
  # - Core CTest
  enable_testing()
  include(CTest)

  # - Geant4_DIR is needed to locate GeantConfig.cmake file required
  # by tests and examples
  set(Geant4_DIR ${PROJECT_BINARY_DIR} CACHE PATH "Current build directory")

  # - Add datasets to testing environment
  geant4_get_datasetnames(_dslist)
  foreach(_ds ${_dslist})
    geant4_get_dataset_property(${_ds} ENVVAR _dsenvvar)
    geant4_get_dataset_property(${_ds} BUILD_DIR _dspath)
    list(APPEND GEANT4_TEST_ENVIRONMENT ${_dsenvvar}=${_dspath})
  endforeach()

  # - Add base URL for test reference files
  set(GEANT4_TEST_REFERENCES_URL "http://cern.ch/geant4-data/stt/references/" CACHE
       STRING "base URL for test reference files")
  mark_as_advanced(GEANT4_TEST_REFERENCES_URL)

  # - Add TOOLS_FONT_PATH if freetype enabled
  if(GEANT4_USE_FREETYPE)
    list(APPEND GEANT4_TEST_ENVIRONMENT TOOLS_FONT_PATH=${PROJECT_SOURCE_DIR}/source/externals/g4tools/fonts)
  endif()

  # - Configure 'G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default)'
  # to use TBB
  if(GEANT4_BUILD_MULTITHREADED AND GEANT4_USE_TBB)
    list(APPEND GEANT4_TEST_ENVIRONMENT G4RUN_MANAGER_TYPE=TBB)
  endif()
endif()

#-----------------------------------------------------------------------
# Add Unit Tests if required
#
if(GEANT4_BUILD_TESTS)
  file(GLOB_RECURSE files RELATIVE ${PROJECT_SOURCE_DIR} source/CMakeLists.txt)
  foreach( file ${files} )
    get_filename_component(path ${file} PATH)
    if(path MATCHES "/test$")
      add_subdirectory(${path})
    endif()
  endforeach()
endif()
