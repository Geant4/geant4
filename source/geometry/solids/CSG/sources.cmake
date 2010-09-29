#------------------------------------------------------------------------------
# sources.cmake
# Module : G4csg
# Package: Geant4.src.G4geometry..G4csg
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.1 2010-09-29 18:42:59 bmorgan Exp $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4csg
    HEADERS
        G4Box.hh
        G4Box.icc
        G4CSGSolid.hh
        G4Cons.hh
        G4Cons.icc
        G4Orb.hh
        G4Orb.icc
        G4Para.hh
        G4Para.icc
        G4Sphere.hh
        G4Sphere.icc
        G4Torus.hh
        G4Torus.icc
        G4Trap.hh
        G4Trap.icc
        G4Trd.hh
        G4Trd.icc
        G4Tubs.hh
        G4Tubs.icc
    SOURCES
        G4Box.cc
        G4CSGSolid.cc
        G4Cons.cc
        G4Orb.cc
        G4Para.cc
        G4Sphere.cc
        G4Torus.cc
        G4Trap.cc
        G4Trd.cc
        G4Tubs.cc
    GRANULAR_DEPENDENCIES
        G4geometrymng
        G4globman
        G4graphics_reps
        G4hepnumerics
        G4intercoms
        G4volumes
    GLOBAL_DEPENDENCIES
        G4global
        G4graphics_reps
        G4intercoms
    LINK_LIBRARIES
)

# List any source specific properties here

