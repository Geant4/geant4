# - Basic setup for testing Geant4 using CMake/CTest
#
# 

#---Do the nexesseray to enable CTest and construct a test environment --------------------------------
if(GEANT4_ENABLE_TESTING)

  enable_testing()
  include(CTest)

  #---Many tests require access to the Geant4 data----------------------------------------------------
  if(GEANT4_INSTALL_DATA)
    set(GEANT4_DATA_DIR ${CMAKE_BINARY_DIR}/data CACHE PATH "Directory where the Geant4 data is located")
  else()
    set(GEANT4_DATA_DIR "" CACHE PATH "Directory where the Geant4 data is located")
  endif()
  
  if(NOT GEANT4_DATA_DIR)
    message(WARNING "GEANT4_DATA_DIR not defined! This may cause many Geant4 tests to fail")
    return()
  endif()

  #---Geant4_DIR is needed to locate GeantConfig.cmake file required by tests and examples-------------
  set(Geant4_DIR ${CMAKE_BINARY_DIR} CACHE PATH "Current build directory")

  #---Define the TEST environment (basically the varibles pointing to DATA files)----------------------
  foreach( tuple "G4LEVELGAMMADATA;PhotonEvaporation"
                 "G4LEDATA;G4EMLOW"
                 "G4RADIOACTIVEDATA;RadioactiveDecay"
                 "G4ELASTICDATA;G4ELASTIC"
                 "NeutronHPCrossSections;G4NDL"
                 "G4NEUTRONHPDATA;G4NDL"
                 "G4ABLADATA;G4ABLA"
                 "G4PIIDATA;G4PII"
                 "G4NEUTRONXSDATA;G4NEUTRONXS")
    list(GET tuple 0 envname)
    list(GET tuple 1 dirname)
    GEANT4_LATEST_VERSION(${GEANT4_DATA_DIR} ${dirname} _result)
    list(APPEND GEANT4_TEST_ENVIRONMENT ${envname}=${_result})
  endforeach()
  
endif()

if(GEANT4_BUILD_TESTS)

  #---Add all directories with Unit Tests---------------------------------------------------------------
  file(GLOB_RECURSE files RELATIVE ${CMAKE_SOURCE_DIR} source/CMakeLists.txt)
  foreach( file ${files} )
    get_filename_component(path ${file} PATH)
    if(path MATCHES "/test$")
      message("PATH = ${path}")
      add_subdirectory(${path})
    endif()
  endforeach()
  
endif()
