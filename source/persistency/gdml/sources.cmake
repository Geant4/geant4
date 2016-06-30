#------------------------------------------------------------------------------
# sources.cmake
# Module : G4gdml
# Package: Geant4.src.G4persistency.G4gdml
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 96673 2016-04-29 12:08:20Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${USOLIDS_INCLUDE_DIRS})

# Need XercesC
include_directories(${XERCESC_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/digits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/event/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/divisions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/Boolean/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/cuts/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/run/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/tracking/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/detector/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4gdml
    HEADERS
        G4GDMLAuxStructType.hh
        G4GDMLEvaluator.hh
        G4GDMLMessenger.hh
        G4GDMLParameterisation.hh
        G4GDMLParser.hh
        G4GDMLParser.icc
        G4GDMLRead.hh
        G4GDMLReadDefine.hh
        G4GDMLReadMaterials.hh
        G4GDMLReadParamvol.hh
        G4GDMLReadSetup.hh
        G4GDMLReadSolids.hh
        G4GDMLReadStructure.hh
        G4GDMLWrite.hh
        G4GDMLWriteDefine.hh
        G4GDMLWriteMaterials.hh
        G4GDMLWriteParamvol.hh
        G4GDMLWriteSetup.hh
        G4GDMLWriteSolids.hh
        G4GDMLWriteStructure.hh
        G4STRead.hh
    SOURCES
        G4GDMLEvaluator.cc
        G4GDMLMessenger.cc
        G4GDMLParameterisation.cc
        G4GDMLParser.cc
        G4GDMLRead.cc
        G4GDMLReadDefine.cc
        G4GDMLReadMaterials.cc
        G4GDMLReadParamvol.cc
        G4GDMLReadSetup.cc
        G4GDMLReadSolids.cc
        G4GDMLReadStructure.cc
        G4GDMLWrite.cc
        G4GDMLWriteDefine.cc
        G4GDMLWriteMaterials.cc
        G4GDMLWriteParamvol.cc
        G4GDMLWriteSetup.cc
        G4GDMLWriteSolids.cc
        G4GDMLWriteStructure.cc
        G4STRead.cc
    GRANULAR_DEPENDENCIES
        G4csg
        G4digits
        G4event
        G4geomBoolean
        G4geomdivision
        G4geometrymng
        G4globman
        G4graphics_reps
        G4hepnumerics
        G4hits
        G4intercoms
        G4materials
        G4navigation
        G4partman
        G4baryons
        G4bosons
        G4leptons
        G4procman
        G4cuts
        G4emutils
        G4run
        G4specsolids
        G4track
        G4tracking
        G4volumes
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4event
        G4geometry
        G4global
        G4graphics_reps
        G4intercoms
        G4materials
        G4particles
        G4processes
        G4run
        G4track
        G4tracking
    LINK_LIBRARIES
        ${XERCESC_LIBRARIES}
)

# List any source specific properties here

