#------------------------------------------------------------------------------
# sources.cmake
# Module : G4phonon
# Package: Geant4.src.G4processes.G4phonon
#
# Sources description for a library.
# Lists the sources and headers of the code explicitly.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/10/2013
#
#
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
GEANT4_DEFINE_MODULE(NAME G4solidstate_channeling
    HEADERS
	G4ChannelingTrackData.hh
        G4ChannelingOptrMultiParticleChangeCrossSection.hh
        G4ChannelingOptrChangeCrossSection.hh
        G4ChannelingMaterialData.hh
        G4Channeling.hh
        G4ChannelingECHARM.hh
    SOURCES
        G4ChannelingTrackData.cc
        G4ChannelingOptrMultiParticleChangeCrossSection.cc
        G4ChannelingOptrChangeCrossSection.cc
        G4ChannelingMaterialData.cc
        G4Channeling.cc
        G4ChannelingECHARM.cc
    GRANULAR_DEPENDENCIES
        G4bosons
        G4baryons
        G4geometrymng
        G4globman
        G4materials
        G4partman
        G4emutils
        G4procman
        G4track
        G4volumes
        G4biasing_mgt
        G4biasing_gen
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

