# - Define useful macros for building and installing Geant4 library targets
# - DEPRECATED: Functionality merged into G4DeveloperAPI/_OLD
#               Retained as many modules still include it
#               Simply forward to new module
# - Include Guard
if(__GEANT4MACROLIBRARYTARGETS_INCLUDED)
  return()
endif()
set(__GEANT4MACROLIBRARYTARGETS_INCLUDED TRUE)

include(G4DeveloperAPI)

