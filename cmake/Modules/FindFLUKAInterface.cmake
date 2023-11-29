# - Check that FLUKA_INTERFACE environment was properly set.
# Note: This file can also be placed in geant4_application/cmake/Modules/FindFLUKAInterface.cmake

message(STATUS "FLUKAInterface_INCLUDE_DIR = $ENV{FLUKAInterface_INCLUDE_DIR}")
message(STATUS "FLUKAInterface_LIBRARIES = $ENV{FLUKAInterface_LIBRARIES}")
if (NOT DEFINED ENV{FLUKAInterface_INCLUDE_DIR})
  message(FATAL_ERROR "\n $FLUKAInterface_INCLUDE_DIR is not set! \n Did you source env_FLUKA_G4_interface.sh? \n")
elseif (NOT DEFINED ENV{FLUKAInterface_LIBRARIES})
  message(FATAL_ERROR "\n $FLUKAInterface_LIBRARIES is not set! \n Did you source env_FLUKA_G4_interface.sh? \n")
endif()


string(REPLACE " " ";" FLUKAInterface_INCLUDE_DIR $ENV{FLUKAInterface_INCLUDE_DIR})
string(REPLACE " " ";" FLUKAInterface_LIBRARIES $ENV{FLUKAInterface_LIBRARIES})


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FLUKAInterface DEFAULT_MSG FLUKAInterface_INCLUDE_DIR FLUKAInterface_LIBRARIES)
