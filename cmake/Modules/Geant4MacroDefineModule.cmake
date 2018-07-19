# - Macros for organizing and specializing code in Geant4 modules
# - DEPRECATED: Functionality merged into G4DeveloperAPI/_OLD
#               Retained as many modules still include it
#               Simply forward to new module
# - Include guard
if(__GEANT4MACRODEFINEMODULE_INCLUDED)
  return()
endif()
set(__GEANT4MACRODEFINEMODULE_INCLUDED YES)

include(G4DeveloperAPI)

