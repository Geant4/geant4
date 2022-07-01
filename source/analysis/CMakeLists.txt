# - G4analysis category build
if(GEANT4_USE_FREETYPE)
  add_definitions(-DTOOLS_USE_FREETYPE)
endif()

include(accumulables/sources.cmake)
include(csv/sources.cmake)
include(factory/sources.cmake)
include(hntools/sources.cmake)
include(management/sources.cmake)
include(root/sources.cmake)
include(xml/sources.cmake)

geant4_add_category(G4analysis MODULES G4accumulables G4csv G4analysisfac G4hntools G4analysismng G4root G4xml)

if(GEANT4_USE_HDF5)
  include(hdf5/sources.cmake)
  geant4_category_modules(G4analysis G4hdf5)
endif()

