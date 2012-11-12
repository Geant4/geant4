# - Basic setup for testing Geant4 using CMake/CTest
#

#-----------------------------------------------------------------------
# Configure CTest and relevant Geant4 settings, if required
#
if(GEANT4_ENABLE_TESTING)
  # - Core CTest
  enable_testing()
  include(CTest)

  # - Geant4_DIR is needed to locate GeantConfig.cmake file required 
  # by tests and examples
  set(Geant4_DIR ${CMAKE_BINARY_DIR} CACHE PATH "Current build directory")

  #-----------------------------------------------------------------------
  # OLD STYLE DATA CONFIGURATION
  # - Configure data location, as most tests will require access to these
  #if(GEANT4_INSTALL_DATA)
  #  set(GEANT4_DATA_DIR ${CMAKE_BINARY_DIR}/data CACHE PATH "Directory where the Geant4 data is located")
  #elseif(NOT "$ENV{GEANT4_DATA_DIR}" STREQUAL "")
  #  set(GEANT4_DATA_DIR  "$ENV{GEANT4_DATA_DIR}" CACHE PATH "Directory where the Geant4 data is located" FORCE)
  #endif()
  
  #if(NOT GEANT4_DATA_DIR)
  #  message(STATUS  "  GEANT4_DATA_DIR not defined! This may cause many Geant4 tests to fail\n"
  #                 "     Add -DGEANT4_DATA_DIR=<path> to the cmake command or define environment")
                  #  return()
                  #endif()

  # - Configure test environment (basically, the data library pointers)
  #foreach( tuple "G4NEUTRONHPDATA;G4NDL"
  #               "G4LEDATA;G4EMLOW"
  #               "G4LEVELGAMMADATA;PhotonEvaporation"
  #               "G4RADIOACTIVEDATA;RadioactiveDecay"
  #               "G4NEUTRONXSDATA;G4NEUTRONXS"
  #               "G4PIIDATA;G4PII"
  #               "G4REALSURFACEDATA;RealSurface"
  #               "G4SAIDXSDATA;G4SAIDDATA"                 
  #               )
  #  list(GET tuple 0 envname)
  #  list(GET tuple 1 dirname)
  #  GEANT4_LATEST_VERSION(${GEANT4_DATA_DIR} ${dirname} _result)
  #  #list(APPEND GEANT4_TEST_ENVIRONMENT ${envname}=${_result})
  #endforeach()

  #-----------------------------------------------------------------------
  # - Configure data using new style API methods
  geant4_get_datasetnames(_dslist)
  foreach(_ds ${_dslist})
    geant4_get_dataset_property(${_ds} ENVVAR _dsenvvar)
    geant4_get_dataset_property(${_ds} BUILD_DIR _dspath)
    list(APPEND GEANT4_TEST_ENVIRONMENT ${_dsenvvar}=${_dspath})
  endforeach()
endif()

#-----------------------------------------------------------------------
# Add Unit Tests if required
#
if(GEANT4_BUILD_TESTS)
  file(GLOB_RECURSE files RELATIVE ${CMAKE_SOURCE_DIR} source/CMakeLists.txt)
  foreach( file ${files} )
    get_filename_component(path ${file} PATH)
    if(path MATCHES "/test$")
      add_subdirectory(${path})
    endif()
  endforeach()
endif()
