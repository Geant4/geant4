#------------------------------------------------------------------------------
# sources.cmake
# Module : G4scoring
# Package: Geant4.src.G4processes.G4scoring
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
GEANT4_DEFINE_MODULE(NAME G4scoring
    HEADERS
        G4EnergySplitter.hh
        G4EnergySplitter.icc
        G4ParallelWorldProcess.hh
        G4ParallelWorldProcessStore.hh
        G4ParallelWorldScoringProcess.hh
        G4ScoreSplittingProcess.hh
    SOURCES
        G4EnergySplitter.cc
        G4ParallelWorldProcess.cc
        G4ParallelWorldProcessStore.cc
        G4ParallelWorldScoringProcess.cc
        G4ScoreSplittingProcess.cc
    GRANULAR_DEPENDENCIES
        G4cuts
        G4detector
        G4muons
        G4emutils
        G4geombias
        G4geometrymng
        G4globman
        G4hits
        G4intercoms
        G4magneticfield
        G4materials
        G4navigation
        G4partman
        G4procman
        G4track
        G4transportation
        G4volumes
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4geometry
        G4global
        G4intercoms
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

