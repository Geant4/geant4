# - Basic setup for testing Geant4 using CTest
# Still rather rough.
# 


#---Deduce the build name--------------------------------------------------------
set(BUILDNAME ${GEANT4_SYSTEM}-${GEANT4_COMPILER}-${CMAKE_BUILD_TYPE})
enable_testing()
include(CTest)


#---Many tests require access to the Geant4 data----------------------------------------------------
if(GEANT4_INSTALL_DATA)
  set(GEANT4_DATA_DIR ${CMAKE_BINARY_DIR}/data CACHE PATH "Directory where the Geant4 data is located")
else()
  set(GEANT4_DATA_DIR "" CACHE PATH "Directory where the Geant4 data is located")
endif()

set(Geant4_DIR ${CMAKE_BINARY_DIR} CACHE PATH "Current build directory")


if(NOT GEANT4_DATA_DIR)
  message(WARNING "GEANT4_DATA_DIR not defined! This may cause many Geant4 tests to fail")
  return()
endif()

#---Helper function to locate the most recent version of a dataset-----------------------------------
function(GEANT4_LATEST_VERSION dir name var)
  file(GLOB files RELATIVE ${dir} ${dir}/${name}*)
  set(newer 0.0)
  foreach(file ${files})
    string(REPLACE ${name} "" version ${file})
    if(${version} VERSION_GREATER ${newer})
      set(newer ${version})
    endif()
  endforeach()
  set(${var} ${dir}/${name}${newer} PARENT_SCOPE)
endfunction()


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
