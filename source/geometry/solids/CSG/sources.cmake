#------------------------------------------------------------------------------
# sources.cmake
# Module : G4csg
# Package: Geant4.src.G4geometry.G4csg
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 103471 2017-04-11 07:31:42Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${USOLIDS_INCLUDE_DIRS})

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
        G4CutTubs.hh
        G4CutTubs.icc
        G4Orb.hh
        G4Orb.icc
        G4OTubs.hh
        G4OTubs.icc
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
        G4UBox.hh
        G4UCons.hh
        G4UOrb.hh
        G4USphere.hh
        G4UTorus.hh
        G4UTrap.hh
        G4UTrd.hh
        G4UTubs.hh
    SOURCES
        G4Box.cc
        G4CSGSolid.cc
        G4Cons.cc
        G4CutTubs.cc
        G4Orb.cc
        G4OTubs.cc
        G4Para.cc
        G4Sphere.cc
        G4Torus.cc
        G4Trap.cc
        G4Trd.cc
        G4Tubs.cc
        G4UBox.cc
        G4UCons.cc
        G4UOrb.cc
        G4USphere.cc
        G4UTorus.cc
        G4UTrap.cc
        G4UTrd.cc
        G4UTubs.cc
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
        ${USOLIDS_LIBRARIES}
)

# List any source specific properties here

