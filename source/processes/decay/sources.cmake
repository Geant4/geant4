#------------------------------------------------------------------------------
# sources.cmake
# Module : G4decay
# Package: Geant4.src.G4processes.G4decay
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
GEANT4_DEFINE_MODULE(NAME G4decay
    HEADERS
        G4Decay.hh
        G4DecayProcessType.hh
        G4DecayWithSpin.hh
        G4PionDecayMakeSpin.hh
        G4UnknownDecay.hh
        G4VExtDecayer.hh
    SOURCES
        G4Decay.cc
        G4DecayWithSpin.cc
        G4PionDecayMakeSpin.cc
        G4UnknownDecay.cc
    GRANULAR_DEPENDENCIES
        G4csg
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

