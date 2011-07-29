#---Many tests require access to the Geant4 data----------------------------------------------------
if(GEANT4_INSTALL_DATA)
  set(GEANT4_DATA_DIR ${CMAKE_BINARY_DIR}/data CACHE PATH "Directory where the Geant4 data is located")
else()
  set(GEANT4_DATA_DIR "" CACHE PATH "Directory where the Geant4 data is located")
endif()

if(NOT GEANT4_DATA_DIR)
  message("GEANT4_DATA_DIR not defined! This may cause many Geant4 tests to fail")
endif()

set(Geant4_DIR ${CMAKE_BINARY_DIR} CACHE PATH "Current build directory")

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
               "G4ABLADATA;G4ABLA")
  list(GET tuple 0 envname)
  list(GET tuple 1 dirname)
  GEANT4_LATEST_VERSION(${GEANT4_DATA_DIR} ${dirname} _result)
  list(APPEND GEANT4_TEST_ENVIRONMENT ${envname}=${_result}) 
endforeach()

#---Add all subdirectories that contains a CMakeLists.txt file -------------------------------------
file(GLOB_RECURSE files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CMakeLists.txt)
foreach( file ${files} )
  get_filename_component(path ${file} PATH)
  if(path)
    add_subdirectory(${path})
  endif()
endforeach()
