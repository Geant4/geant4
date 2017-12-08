#------------------------------------------------------------------------------
# sources.cmake
# Module : G4ions
# Package: Geant4.src.G4particles..G4ions
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 106143 2017-09-14 06:34:42Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4ions
    HEADERS
        G4Alpha.hh
        G4AntiAlpha.hh
        G4AntiDeuteron.hh
        G4AntiHe3.hh
        G4AntiTriton.hh
        G4Deuteron.hh
        G4GenericIon.hh
        G4He3.hh
        G4IonConstructor.hh
        G4Triton.hh
        G4GenericMuonicAtom.hh
    SOURCES
        G4Alpha.cc
        G4AntiAlpha.cc
        G4AntiDeuteron.cc
        G4AntiHe3.cc
        G4AntiTriton.cc
        G4Deuteron.cc
        G4GenericIon.cc
        G4He3.cc
        G4IonConstructor.cc
        G4Triton.cc
        G4GenericMuonicAtom.cc
    GRANULAR_DEPENDENCIES
        G4globman
        G4materials
        G4partman
    GLOBAL_DEPENDENCIES
        G4global
        G4materials
    LINK_LIBRARIES
)

# List any source specific properties here

