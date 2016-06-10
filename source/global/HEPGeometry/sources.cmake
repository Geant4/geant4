#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hepgeometry
# Package: Geant4.src.G4global.G4hepgeometry
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 66892 2013-01-17 10:57:59Z gunter $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hepgeometry
    HEADERS
        G4LorentzRotation.hh
        G4LorentzVector.hh
        G4Normal3D.hh
        G4Plane3D.hh
        G4Point3D.hh
        G4Transform3D.hh
        G4Vector3D.hh
        geomdefs.hh
    SOURCES
    GRANULAR_DEPENDENCIES
        G4globman
    GLOBAL_DEPENDENCIES
    LINK_LIBRARIES
)

# List any source specific properties here

