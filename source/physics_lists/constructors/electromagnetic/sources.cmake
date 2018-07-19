#------------------------------------------------------------------------------
# sources.cmake
# Module : G4phys_ctor_em
# Package: Geant4.src.G4physicslists.G4physlist_ctors.G4physlist_ctor_em
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 10/01/2013
#
# $Id: sources.cmake 104020 2017-05-08 07:34:58Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/digits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/event/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/magneticfield/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/shortlived/include)
include_directories(${CMAKE_SOURCE_DIR}/source/physics_lists/builders/include)
include_directories(${CMAKE_SOURCE_DIR}/source/physics_lists/constructors/factory/include)
include_directories(${CMAKE_SOURCE_DIR}/source/physics_lists/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/cuts/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/decay/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/processes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/models/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/molecules/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/molecules/types/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/highenergy/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/lowenergy/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/muons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/standard/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/xrays/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/optical/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/transportation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/run/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/tracking/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4phys_ctor_em
  HEADERS
    G4EmDNAChemistry.hh
    G4EmDNAPhysics.hh
    G4EmDNAPhysics_option1.hh
    G4EmDNAPhysics_option2.hh
    G4EmDNAPhysics_option3.hh
    G4EmDNAPhysics_option4.hh
    G4EmDNAPhysics_option5.hh
    G4EmDNAPhysics_option6.hh
    G4EmDNAPhysics_option7.hh
    G4EmDNAPhysics_stationary.hh
    G4EmDNAPhysics_stationary_option2.hh
    G4EmDNAPhysics_stationary_option4.hh
    G4EmDNAPhysics_stationary_option6.hh
    G4EmDNAPhysicsActivator.hh
    G4EmLEPTSPhysics.hh
    G4EmLivermorePhysics.hh
    G4EmLivermorePolarizedPhysics.hh
    G4EmLowEPPhysics.hh
    G4EmModelActivator.hh
    G4EmParticleList.hh
    G4EmPenelopePhysics.hh
    G4EmStandardPhysics.hh
    G4EmStandardPhysicsGS.hh
    G4EmStandardPhysicsSS.hh
    G4EmStandardPhysicsWVI.hh
    G4EmStandardPhysics_option1.hh
    G4EmStandardPhysics_option2.hh
    G4EmStandardPhysics_option3.hh
    G4EmStandardPhysics_option4.hh
    G4OpticalPhysics.hh
    G4OpticalPhysicsMessenger.hh
    G4OpticalProcessIndex.hh
  SOURCES
    G4EmDNAChemistry.cc
    G4EmDNAPhysics.cc
    G4EmDNAPhysics_option1.cc
    G4EmDNAPhysics_option2.cc
    G4EmDNAPhysics_option3.cc
    G4EmDNAPhysics_option4.cc
    G4EmDNAPhysics_option5.cc
    G4EmDNAPhysics_option6.cc
    G4EmDNAPhysics_option7.cc
    G4EmDNAPhysics_stationary.cc
    G4EmDNAPhysics_stationary_option2.cc
    G4EmDNAPhysics_stationary_option4.cc
    G4EmDNAPhysics_stationary_option6.cc
    G4EmDNAPhysicsActivator.cc
    G4EmLEPTSPhysics.cc
    G4EmLivermorePhysics.cc
    G4EmLivermorePolarizedPhysics.cc
    G4EmLowEPPhysics.cc
    G4EmModelActivator.cc
    G4EmParticleList.cc
    G4EmPenelopePhysics.cc
    G4EmStandardPhysics.cc
    G4EmStandardPhysicsGS.cc
    G4EmStandardPhysicsSS.cc
    G4EmStandardPhysicsWVI.cc
    G4EmStandardPhysics_option1.cc
    G4EmStandardPhysics_option2.cc
    G4EmStandardPhysics_option3.cc
    G4EmStandardPhysics_option4.cc
    G4OpticalPhysics.cc
    G4OpticalPhysicsMessenger.cc
  GRANULAR_DEPENDENCIES
    G4baryons
    G4bosons
    G4cuts
    G4decay
    G4digits
    G4emdna-utils
    G4emdna-processes
    G4emdna-molman
    G4emdna-moltypes
    G4emdna-models
    G4emhighenergy
    G4emlowenergy
    G4emstandard
    G4emutils
    G4event
    G4geometrymng
    G4globman
    G4hits
    G4intercoms
    G4ions
    G4leptons
    G4magneticfield
    G4materials
    G4mesons
    G4muons
    G4navigation
    G4optical
    G4partman
    G4phys_builders
    G4phys_ctor_factory
    G4physlist_util
    G4procman
    G4run
    G4shortlived
    G4track
    G4tracking
    G4transportation
    G4volumes
    G4xrays
  GLOBAL_DEPENDENCIES
    G4digits_hits
    G4event
    G4geometry
    G4global
    G4intercoms
    G4materials
    G4particles
    G4processes
    G4run
    G4track
    G4tracking
  LINK_LIBRARIES
)

# List any source specific properties here

