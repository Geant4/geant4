#---Adding all advanced examples subdirectories explicitly 

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

#----------------------------------------------------------------------------
# some examples require Geant4 build with optional packages
#
find_package(Geant4)

find_package(ROOT QUIET)

add_subdirectory(air_shower)
add_subdirectory(amsEcal)
add_subdirectory(brachytherapy)

if(Geant4_gdml_FOUND)
  add_subdirectory(ChargeExchangeMC)
else()  
  message(STATUS "G4 Examples: GDML not installed --> ChargeExchangeMC example disabled")  
endif()

if(NOT Geant4_builtin_clhep_FOUND)
  add_subdirectory(nanobeam)
else()
  message(STATUS "G4 Examples: CLHEP external package not found  --> nanobeam disabled")
endif() 

add_subdirectory(composite_calorimeter)

if(ROOT_FOUND)
  add_subdirectory(doiPET)
else()
  message(STATUS "G4 Examples: ROOT package not found --> doiPET disabled")
endif()  

add_subdirectory(eRosita)
add_subdirectory(gammaknife)
add_subdirectory(gammaray_telescope)
add_subdirectory(hadrontherapy)
add_subdirectory(human_phantom)
add_subdirectory(iort_therapy)
add_subdirectory(lAr_calorimeter)
add_subdirectory(medical_linac)
add_subdirectory(microbeam)
add_subdirectory(microelectronics)
add_subdirectory(purging_magnet)
add_subdirectory(underground_physics)
add_subdirectory(xray_fluorescence)
add_subdirectory(STCyclotron)

