#------------------------------------------------------------------------------
# sources.cmake
# Module : G4partman
# Package: Geant4.src.G4particles.G4partman
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
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4partman
    HEADERS
        G4DalitzDecayChannel.hh
        G4DecayProducts.hh
        G4DecayTable.hh
        G4DecayTableMessenger.hh
        G4DynamicParticle.hh
        G4DynamicParticle.icc
        G4DynamicParticleFastVector.hh
        G4ElectronOccupancy.hh
        G4HyperNucleiProperties.hh
        G4IonTable.hh
        G4Ions.hh
        G4IsotopeProperty.hh
        G4KL3DecayChannel.hh
        G4MuonicAtom.hh
        G4MuonicAtomHelper.hh
        G4MuonDecayChannel.hh
        G4MuonDecayChannelWithSpin.hh
        G4MuonRadiativeDecayChannelWithSpin.hh
        G4NeutronBetaDecayChannel.hh
        G4NucleiProperties.hh
        G4NucleiPropertiesTableAME03.hh
        G4NucleiPropertiesTableAME12.hh
        G4NucleiPropertiesTheoreticalTable.hh
        G4NuclideTable.hh
        G4NuclideTableMessenger.hh
        G4PDGCodeChecker.hh
        G4PDefManager.hh
        G4ParticleDefinition.hh
        G4ParticleDefinition.icc
        G4ParticleMessenger.hh
        G4ParticleMomentum.hh
        G4ParticlePropertyData.hh
        G4ParticlePropertyData.icc
        G4ParticlePropertyMessenger.hh
        G4ParticlePropertyTable.hh
        G4ParticleTable.hh
        G4ParticleTable.icc
        G4ParticleTableIterator.hh
        G4ParticleWithCuts.hh
        G4ParticlesWorkspace.hh
        G4PhaseSpaceDecayChannel.hh
        G4PionRadiativeDecayChannel.hh
        G4PrimaryParticle.hh
        G4PrimaryVertex.hh
        G4TauLeptonicDecayChannel.hh
        G4VDecayChannel.hh
        G4VIsotopeTable.hh
        G4VUserPrimaryParticleInformation.hh
        G4VUserPrimaryVertexInformation.hh
        pwdefs.hh
	SOURCES
        G4DalitzDecayChannel.cc
        G4DecayProducts.cc
        G4DecayTable.cc
        G4DecayTableMessenger.cc
        G4DynamicParticle.cc
        G4ElectronOccupancy.cc
        G4HyperNucleiProperties.cc
        G4IonTable.cc
        G4Ions.cc
        G4IsotopeProperty.cc
        G4KL3DecayChannel.cc
        G4MuonicAtom.cc
        G4MuonicAtomHelper.cc
        G4MuonDecayChannel.cc
        G4MuonDecayChannelWithSpin.cc
        G4MuonRadiativeDecayChannelWithSpin.cc
        G4NeutronBetaDecayChannel.cc
        G4NucleiProperties.cc
        G4NucleiPropertiesTableAME03.cc
        G4NucleiPropertiesTableAME12.cc
        G4NucleiPropertiesTheoreticalTableA.cc
        G4NucleiPropertiesTheoreticalTableB.cc
        G4NuclideTable.cc
        G4NuclideTableMessenger.cc
        G4PDGCodeChecker.cc
        G4PDefManager.cc
        G4ParticleDefinition.cc
        G4ParticleMessenger.cc
        G4ParticlePropertyData.cc
        G4ParticlePropertyMessenger.cc
        G4ParticlePropertyTable.cc
        G4ParticleTable.cc
        G4ParticlesWorkspace.cc
        G4PhaseSpaceDecayChannel.cc
        G4PionRadiativeDecayChannel.cc
        G4PrimaryParticle.cc
        G4PrimaryVertex.cc
        G4TauLeptonicDecayChannel.cc
        G4VDecayChannel.cc
        G4VIsotopeTable.cc
        G4VUserPrimaryParticleInformation.cc
        G4VUserPrimaryVertexInformation.cc
	GRANULAR_DEPENDENCIES
        G4geometrymng
        G4globman
        G4intercoms
        G4materials
	GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4intercoms
        G4materials
	LINK_LIBRARIES
	)

      # List any source specific properties here
      if("$ENV{G4NucleiProperties_USE_OLD_AME_TABLE}")
	set_source_files_properties(
	  ${G4partman_SOURCES}
	  PROPERTIES COMPILE_DEFINITIONS G4NucleiProperties_USE_OLD_AME_TABLE
	  )
      endif()


