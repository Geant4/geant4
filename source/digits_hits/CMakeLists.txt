# - G4digits_hits category build
add_definitions(-DG4DIGI_ALLOC_EXPORT)

include(detector/sources.cmake)
include(digits/sources.cmake)
include(hits/sources.cmake)
include(scorer/sources.cmake)
include(utils/sources.cmake)

geant4_add_category(G4digits_hits MODULES G4detector G4digits G4hits G4detscorer G4detutils)
