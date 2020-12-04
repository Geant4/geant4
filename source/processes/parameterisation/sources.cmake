#------------------------------------------------------------------------------
# sources.cmake
# Module : G4parameterisation
# Package: Geant4.src.G4processes.G4parameterisation
#
# Sources description for a library.
# Lists the sources and headers of the code explicitly.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
#
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
GEANT4_DEFINE_MODULE(NAME G4parameterisation
    HEADERS
        G4FastHit.hh
        G4FastSimHitMaker.hh
        G4FastSimulationHelper.hh
        G4FastSimulationManager.hh
        G4FastSimulationManagerProcess.hh
        G4FastSimulationMessenger.hh
        G4FastSimulationProcessType.hh
        G4FastSimulationVector.hh
        G4FastSimulationVector.icc
        G4FastStep.hh
        G4FastStep.icc
        G4FastTrack.hh
        G4GlobalFastSimulationManager.hh
        G4VFastSimSensitiveDetector.hh
        G4VFastSimulationModel.hh
    SOURCES
        G4FastHit.cc
        G4FastSimHitMaker.cc
        G4FastSimulationHelper.cc
        G4FastSimulationManager.cc
        G4FastSimulationManagerProcess.cc
        G4FastSimulationMessenger.cc
        G4FastStep.cc
        G4FastTrack.cc
        G4GlobalFastSimulationManager.cc
        G4VFastSimulationModel.cc
    GRANULAR_DEPENDENCIES
        G4geometrymng
        G4globman
        G4intercoms
        G4magneticfield
        G4materials
        G4navigation
        G4partman
        G4procman
        G4track
        G4volumes
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4intercoms
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

